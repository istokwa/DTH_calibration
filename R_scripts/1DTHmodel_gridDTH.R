# ============================================================
# FULL SCRIPT (keeps your Daylength + Temperature blocks as-is)
# Goal: Compute predicted DTH for ALL parameter combinations in bigspace
# Output: bigspace_with_DTH.xlsx with MODEL1..MODEL<number_of_sites>
# ============================================================

library(Rcpp)
library(foreach)
library(doParallel)
library(readxl)
library(writexl)
library(tictoc)
library(RcppParallel)
library(suncalc)
library(solartime)

setwd("D:/NU/Repository/GP_DTH-Rice/")

# Use your compiled C++ (must export calculateDTH and/or calculateDTH_optimized)
sourceCpp("src_cpp/mse_v05.cpp")

phenotypes <- read_excel("data/raw/phenotypes_DVIDVRv02.xlsx", sheet = "Sheet2")

experimentID <- "(01-14-2026)-finer-002"

# ============================================================
# Generate DTH (observed)
# ============================================================
dth <- data.frame(phenotypes[, 1])
colnames(dth) <- "taxa"

number_of_sites <- (ncol(phenotypes) - 1) / 3

for (i in 1:number_of_sites) {
  oldcolnames <- colnames(dth)
  
  start_date <- as.Date(phenotypes[[i * 3]])
  end_date   <- as.Date(phenotypes[[i * 3 + 1]])
  
  dth[, i + 1] <- as.numeric(
    difftime(end_date, start_date, units = "days")
  )
  
  colnames(dth) <- c(oldcolnames, paste("dth", i))
}

# ============================================================
# Generate Daylengths  (YOUR CODE — UNCHANGED)
# ============================================================
daylegnths <- matrix(data = 0, nrow = 150, ncol = number_of_sites)

for (i in 1:number_of_sites) {
  if (phenotypes[1, i * 3 - 1] == "ishigaki") {
    latitude <- 24.3784258  # JIRCAS Ishigaki
    longitude <- 124.1952777
    elevation <- 32
  }
  
  if (phenotypes[1, i * 3 - 1] == "nagoya") {
    latitude <- 35.11167  # Togo
    longitude <- 137.08195
    elevation <- 65
  }
  
  if (phenotypes[1, i * 3 - 1] == "iwate") {
    latitude <- 39.35107
    longitude <- 141.10688
    elevation <- 63
  }
  
  day0 <- as.Date(phenotypes[[1, i * 3]])
  
  for (j in 1:150) {
    date <- as.Date(day0 + j - 2) # sowing date of the 1st taxa + j
    
    sun_times <- getSunlightTimes(
      date = date,
      lat = latitude,
      lon = longitude,
      keep = c("sunrise", "sunset")
    )
    
    horizon_correction <- -2.076 * sqrt(elevation) / 60
    
    sunrise <- as.POSIXct(sun_times$sunrise, format = "%Y-%m-%dT%H:%M:%S", tz = "UTC")
    sunset <- as.POSIXct(sun_times$sunset, format = "%Y-%m-%dT%H:%M:%S", tz = "UTC")
    
    sunrise_corrected <- sunrise - horizon_correction * 60
    sunset_corrected <- sunset + horizon_correction * 60
    
    daylegnths[j, i] <- sunset_corrected - sunrise_corrected
  }
}

daylengths <- data.frame(daylegnths)

# ============================================================
# Generate Temperature  (YOUR CODE — UNCHANGED)
# ============================================================
temperatures <- NULL

for (i in 1:number_of_sites) {
  file_name <- paste0(phenotypes[1, i * 3 - 1], ".xlsx")
  file_path <- file.path("data", "raw", file_name)
  
  temperaturedata <- read_excel(file_path)
  
  day0 <- phenotypes[[1, i * 3]]
  row_of_day0 <- which(temperaturedata$Date == as.Date(day0))
  temp <- temperaturedata[row_of_day0:(row_of_day0 + 150), 2]
  
  temperatures <- if (is.null(temperatures)) temp else cbind(temperatures, temp)
}

# ============================================================
# IMPORTANT: Convert daylengths to NUMERIC HOURS
# Your current daylengths are difftime objects (seconds)
# ============================================================
temperatures <- as.data.frame(temperatures)

# Optional: name the columns
colnames(daylengths)   <- paste0("DL", 1:number_of_sites)
colnames(temperatures) <- paste0("T",  1:number_of_sites)

# ============================================================
# DTH calculator for one parameter set (returns vector across sites)
# ============================================================
calcDTH <- function(temperatures,
                    daylengths,
                    G,
                    Th,
                    Lc,
                    A,
                    B,
                    DVIstar) {
  
  n_sites <- ncol(temperatures)
  out <- numeric(n_sites)
  
  use_fast <- exists("calculateDTH_optimized")
  
  for (i in 1:n_sites) {
    if (use_fast) {
      out[i] <- calculateDTH_optimized(
        as.numeric(temperatures[, i]),
        as.numeric(daylengths[, i]),
        G, Th, Lc, A, B, DVIstar
      )
    } else {
      out[i] <- calculateDTH(
        as.numeric(temperatures[, i]),
        as.numeric(daylengths[, i]),
        G, Th, Lc, A, B, DVIstar
      )
    }
  }
  
  # ✅ force whole-number day counts
  out <- as.integer(round(out))  # use round(); or ceiling(); or floor()
  
  return(out)
}

# ============================================================
# Coarse Space (your parameter ranges)
# ============================================================
bigspace <- list()
counter <- 1

for (G in seq(50, 60, by = 10)) {
  for (Th in seq(15, 20, by = 5)){
    for (Lc in seq(15, 15, by = 0)){
      for (A in seq(0, 1, by = 0.01)){
        for (B in seq(7, 15, by = 0.02)) {
          for (DVIstar in seq(0.4, 0.8, by = 0.05)) {
            bigspace[[counter]] <- c(G, Th, Lc, A, B, DVIstar)
            counter <- counter + 1
          }
        }
      }
    }
  }
}

bigspace <- as.data.frame(do.call(rbind, bigspace))
colnames(bigspace) <- c("G", "Th", "Lc", "A", "B", "DVIstar")

cat("Total parameter combinations:", nrow(bigspace), "\n")
cat("Number of environments/sites:", number_of_sites, "\n")

# ============================================================
# Predict DTH for ALL combinations (ALL environments)
# ============================================================
tic("Predicting DTH for all combinations...")

pred_mat <- t(apply(bigspace, 1, function(p) {
  calcDTH(
    temperatures = temperatures,
    daylengths   = daylengths,
    G       = as.numeric(p["G"]),
    Th      = as.numeric(p["Th"]),
    Lc      = as.numeric(p["Lc"]),
    A       = as.numeric(p["A"]),
    B       = as.numeric(p["B"]),
    DVIstar = as.numeric(p["DVIstar"])
  )
}))

toc()

colnames(pred_mat) <- paste0("MODEL", 1:number_of_sites)

bigspace_with_DTH <- cbind(bigspace, pred_mat)

# ============================================================
# Save output
# ============================================================
output_dir <- file.path("results", "(04) Bigspace_DTH")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

outfile <- file.path(output_dir, paste0(experimentID,".xlsx"))
write_xlsx(list(AllCombos_PredDTH = as.data.frame(bigspace_with_DTH)), outfile)

cat("✅ Saved:", outfile, "\n")
