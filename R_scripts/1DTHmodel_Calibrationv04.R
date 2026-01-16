# ============================================================
# 1DTHmodel_Calibrationv03.R  (FULL SCRIPT - Coarse + Fine + Finer)
# - Fine scan runs on ALL coarse exact minima seeds
# - Finer scan runs on ALL fine exact minima seeds
# - Predicted DTH is rounded to integer BEFORE computing MSE
# - Saves Coarse/Fine/Finer exact minima to one Excel file
# ============================================================

suppressPackageStartupMessages({
  library(Rcpp)
  library(foreach)
  library(doParallel)
  library(readxl)
  library(writexl)
  library(tictoc)
  library(RcppParallel)
  library(suncalc)
  library(solartime)
})

setwd("D:/NU/Repository/GP_DTH-Rice/")

sourceCpp("src_cpp/mse_v05.cpp")

phenotypes <- read_excel("data/raw/phenotypes_DVIDVRv02.xlsx", sheet = "Sheet3")

experimentID <- "(01-15-2026)-BASE-007"

# ============================================================
# Generate DTH
# ============================================================
dth <- data.frame(phenotypes[, 1])
colnames(dth) <- "taxa"

number_of_sites <- (ncol(phenotypes) - 1) / 3

for (i in 1:number_of_sites) {
  oldcolnames <- colnames(dth)
  
  start_date <- as.Date(phenotypes[[i * 3]])
  end_date   <- as.Date(phenotypes[[i * 3 + 1]])
  
  dth[, i + 1] <- as.numeric(difftime(end_date, start_date, units = "days"))
  
  colnames(dth) <- c(oldcolnames, paste0("dth", i))
}

# ============================================================
# Generate Daylengths
# ============================================================
daylegnths <- matrix(data = 0, nrow = 151, ncol = number_of_sites)

for (i in 1:number_of_sites) {
  if (phenotypes[1, i * 3 - 1] == "ishigaki") {
    latitude <- 24.3784258
    longitude <- 124.1952777
    elevation <- 32
  }
  
  if (phenotypes[1, i * 3 - 1] == "nagoya") {
    latitude <- 35.11167
    longitude <- 137.08195
    elevation <- 65
  }
  
  if (phenotypes[1, i * 3 - 1] == "iwate") {
    latitude <- 39.35107
    longitude <- 141.10688
    elevation <- 64
  }
  
  day0 <- as.Date(phenotypes[[1, i * 3]])
  
  for (j in 1:151) {
    date <- as.Date(day0 + j - 2)
    
    sun_times <- getSunlightTimes(
      date = date,
      lat = latitude,
      lon = longitude,
      keep = c("sunrise", "sunset")
    )
    
    horizon_correction <- -2.076 * sqrt(elevation) / 60
    
    sunrise <- as.POSIXct(sun_times$sunrise, format = "%Y-%m-%dT%H:%M:%S", tz = "UTC")
    sunset  <- as.POSIXct(sun_times$sunset,  format = "%Y-%m-%dT%H:%M:%S", tz = "UTC")
    
    sunrise_corrected <- sunrise - horizon_correction * 60
    sunset_corrected  <- sunset  + horizon_correction * 60
    
    daylegnths[j, i] <- sunset_corrected - sunrise_corrected
  }
}

daylengths <- data.frame(daylegnths)

# ============================================================
# Generate Temperature
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
# Coarse Space
# ============================================================
bigspace <- list()
counter <- 1

for (G in seq(70, 20, by = -15)) {
  for (Th in seq(10, 30, by = 10)) {
    for (Lc in seq(10, 15, by = 1)) {
      for (A in seq(0, 2, by = 0.1)) {
        for (B in seq(0, 15, by = 0.5)) {
          for (DVIstar in seq(0.2, 0.6, by = 0.2)) {
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

# ============================================================
# Fine Space (offset template)
# ============================================================
template_of_smallspace <- list()
counter <- 1

for (G in seq(10, -10, by = -5)) {
  for (Th in seq(0, 0, by = 0)) {
    for (Lc in seq(0, 0, by = 0)) {
      for (A in seq(0.2, -0.2, by = -0.1)) {
        for (B in seq(2, -2, by = -1)) {
          for (DVIstar in seq(0.1, -0.1, by = -0.05)) {
            template_of_smallspace[[counter]] <- c(G, Th, Lc, A, B, DVIstar)
            counter <- counter + 1
          }
        }
      }
    }
  }
}

template_of_smallspace <- as.data.frame(do.call(rbind, template_of_smallspace))
colnames(template_of_smallspace) <- c("G", "Th", "Lc", "A", "B", "DVIstar")

# ============================================================
# Finer Space (offset template)
# ============================================================
template_of_smallerspace <- list()
counter <- 1

for (G in seq(2, -2, by = -0.5)) {
  for (Th in seq(0, 0, by = 0)) {
    for (Lc in seq(0, 0, by = 0)) {
      for (A in seq(0.05, -0.05, by = -0.01)) {
        for (B in seq(0.5, -0.5, by = -0.1)) {
          for (DVIstar in seq(0, 0, by = 0)) {
            template_of_smallerspace[[counter]] <- c(G, Th, Lc, A, B, DVIstar)
            counter <- counter + 1
          }
        }
      }
    }
  }
}

template_of_smallerspace <- as.data.frame(do.call(rbind, template_of_smallerspace))
colnames(template_of_smallerspace) <- c("G", "Th", "Lc", "A", "B", "DVIstar")

# ============================================================
# DF for results
# ============================================================
DVRparams <- data.frame(matrix(ncol = 8 + 2 * number_of_sites, nrow = 0))
params <- vector("list", length = 8 + 2 * number_of_sites)

columnnames <- c("Taxa", "G", "Th", "Lc", "A", "B", "DVIstar", "MSE")
for (i in 1:number_of_sites) columnnames <- c(columnnames, paste0("DTH", i))
for (i in 1:number_of_sites) columnnames <- c(columnnames, paste0("MODEL", i))

colnames(DVRparams) <- columnnames
names(params) <- columnnames

all_exact_minima_coarse <- data.frame()
all_exact_minima_fine   <- data.frame()
all_exact_minima_finer  <- data.frame()

# ============================================================
# Helpers
# ============================================================
add_obs_pred_cols <- function(df, tempdth, model, number_of_sites) {
  tempdth <- as.numeric(tempdth)
  model   <- as.numeric(model)
  
  for (k in 1:number_of_sites) df[[paste0("DTH", k)]]   <- tempdth[k]
  for (k in 1:number_of_sites) df[[paste0("MODEL", k)]] <- model[k]
  df
}

make_best_row <- function(temp_best, min_mse, TaxaName, tempdth, model, number_of_sites) {
  best_row <- data.frame(
    temp_best,
    MSE  = min_mse,
    Taxa = as.character(TaxaName[[1]]),
    stringsAsFactors = FALSE
  )
  add_obs_pred_cols(best_row, tempdth, model, number_of_sites)
}

# ============================================================
# Scan taxa by taxa
# ============================================================
param_cols <- c("G", "Th", "Lc", "A", "B", "DVIstar")

for (taxa in 1:3) {
  tic()
  
  TaxaName <- phenotypes[taxa, 1]
  tempdth <- dth[taxa, 2:(number_of_sites + 1)]
  newMSE <- 99999.0
  
  if (!any(is.na(tempdth))) {
    
    # ----------------------------------------------------------
    # MSE evaluator (INTEGER predicted DTH)
    # ----------------------------------------------------------
    getanMSE <- function(par) {
      G <- par[1]; Th <- par[2]; Lc <- par[3]; A <- par[4]; B <- par[5]; DVIstar <- par[6]
      
      pred <- numeric(number_of_sites)
      for (i in 1:number_of_sites) {
        pred[i] <- calculateDTH_optimized(
          temperatures[, i],
          daylengths[, i],
          G, Th, Lc, A, B, DVIstar
        )
      }
      
      pred <- round(pred)
      mean((as.numeric(tempdth) - pred)^2)
    }
    
    # ----------------------------------------------------------
    # Coarse scan
    # ----------------------------------------------------------
    cat("Coarse scan for", as.character(TaxaName[[1]]), "\n")
    
    results <- as.numeric(apply(bigspace[, param_cols, drop = FALSE], 1, getanMSE))
    min_index <- which.min(results)
    
    temp_best <- bigspace[min_index, , drop = FALSE]
    min_mse <- results[min_index]
    
    results_df <- data.frame(bigspace, MSE = results)
    equal_minima <- results_df[results_df$MSE == min_mse, , drop = FALSE]
    equal_minima$Taxa <- as.character(TaxaName[[1]])
    
    # build MODEL for coarse best (integer)
    model_coarse <- numeric(number_of_sites)
    for (i in 1:number_of_sites) {
      model_coarse[i] <- round(
        calculateDTH_optimized(
          temperatures[, i],
          daylengths[, i],
          temp_best$G, temp_best$Th, temp_best$Lc, temp_best$A, temp_best$B, temp_best$DVIstar
        )
      )
    }
    
    equal_minima <- add_obs_pred_cols(equal_minima, tempdth, model_coarse, number_of_sites)
    all_exact_minima_coarse <- rbind(all_exact_minima_coarse, equal_minima)
    
    # ----------------------------------------------------------
    # Fine scan (ALL coarse exact minima seeds)
    # ----------------------------------------------------------
    cat("Fine scan for", as.character(TaxaName[[1]]), "\n")
    
    seeds <- equal_minima[, param_cols, drop = FALSE]
    
    fine_results_all <- NULL
    best_fine_mse <- Inf
    best_fine_row <- NULL
    
    for (s in 1:nrow(seeds)) {
      seed_vec <- as.numeric(seeds[s, param_cols])
      
      smallspace_seed <- template_of_smallspace
      smallspace_seed[, param_cols] <- sweep(
        template_of_smallspace[, param_cols, drop = FALSE],
        2,
        seed_vec,
        FUN = "+"
      )
      
      results_seed <- as.numeric(apply(smallspace_seed[, param_cols, drop = FALSE], 1, getanMSE))
      results_seed_df <- data.frame(smallspace_seed[, param_cols, drop = FALSE], MSE = results_seed)
      
      results_seed_df$Taxa <- as.character(TaxaName[[1]])
      results_seed_df$SeedID <- s
      
      fine_results_all <- rbind(fine_results_all, results_seed_df)
      
      min_idx_seed <- which.min(results_seed)
      min_mse_seed <- results_seed[min_idx_seed]
      
      if (min_mse_seed < best_fine_mse) {
        best_fine_mse <- min_mse_seed
        best_fine_row <- results_seed_df[min_idx_seed, , drop = FALSE]
      }
    }
    
    equal_minima_fine <- fine_results_all[fine_results_all$MSE == best_fine_mse, , drop = FALSE]
    
    # build MODEL for fine best (integer)
    model_fine <- numeric(number_of_sites)
    for (i in 1:number_of_sites) {
      model_fine[i] <- round(
        calculateDTH_optimized(
          temperatures[, i],
          daylengths[, i],
          best_fine_row$G, best_fine_row$Th, best_fine_row$Lc, best_fine_row$A, best_fine_row$B, best_fine_row$DVIstar
        )
      )
    }
    
    equal_minima_fine <- add_obs_pred_cols(equal_minima_fine, tempdth, model_fine, number_of_sites)
    all_exact_minima_fine <- rbind(all_exact_minima_fine, equal_minima_fine)
    
    # ----------------------------------------------------------
    # Finer scan (ALL fine exact minima seeds)
    # ----------------------------------------------------------
    cat("Finer scan for", as.character(TaxaName[[1]]), "\n")
    
    seeds_fine <- equal_minima_fine[, param_cols, drop = FALSE]
    
    finer_results_all <- NULL
    best_finer_mse <- Inf
    best_finer_row <- NULL
    
    for (s in 1:nrow(seeds_fine)) {
      seed_vec <- as.numeric(seeds_fine[s, param_cols])
      
      finerspace_seed <- template_of_smallerspace
      finerspace_seed[, param_cols] <- sweep(
        template_of_smallerspace[, param_cols, drop = FALSE],
        2,
        seed_vec,
        FUN = "+"
      )
      
      results_seed <- as.numeric(apply(finerspace_seed[, param_cols, drop = FALSE], 1, getanMSE))
      results_seed_df <- data.frame(finerspace_seed[, param_cols, drop = FALSE], MSE = results_seed)
      
      results_seed_df$Taxa <- as.character(TaxaName[[1]])
      results_seed_df$SeedID <- s
      
      finer_results_all <- rbind(finer_results_all, results_seed_df)
      
      min_idx_seed <- which.min(results_seed)
      min_mse_seed <- results_seed[min_idx_seed]
      
      if (min_mse_seed < best_finer_mse) {
        best_finer_mse <- min_mse_seed
        best_finer_row <- results_seed_df[min_idx_seed, , drop = FALSE]
      }
    }
    
    equal_minima_finer <- finer_results_all[finer_results_all$MSE == best_finer_mse, , drop = FALSE]
    
    # build MODEL for finer best (integer)
    model_finer <- numeric(number_of_sites)
    for (i in 1:number_of_sites) {
      model_finer[i] <- round(
        calculateDTH_optimized(
          temperatures[, i],
          daylengths[, i],
          best_finer_row$G, best_finer_row$Th, best_finer_row$Lc, best_finer_row$A, best_finer_row$B, best_finer_row$DVIstar
        )
      )
    }
    
    equal_minima_finer <- add_obs_pred_cols(equal_minima_finer, tempdth, model_finer, number_of_sites)
    all_exact_minima_finer <- rbind(all_exact_minima_finer, equal_minima_finer)
    
    # ----------------------------------------------------------
    # Save DVRparams row for this taxa (use FINER best)
    # ----------------------------------------------------------
    optG       <- best_finer_row$G
    optTh      <- best_finer_row$Th
    optLc      <- best_finer_row$Lc
    optA       <- best_finer_row$A
    optB       <- best_finer_row$B
    optDVIstar <- best_finer_row$DVIstar
    newMSE     <- best_finer_row$MSE
    
    params[["Taxa"]] <- as.character(TaxaName[[1]])
    params[["G"]] <- optG
    params[["Th"]] <- optTh
    params[["Lc"]] <- optLc
    params[["A"]] <- optA
    params[["B"]] <- optB
    params[["DVIstar"]] <- optDVIstar
    params[["MSE"]] <- newMSE
    
    for (i in 1:number_of_sites) params[[paste0("DTH", i)]] <- as.numeric(tempdth[i])
    for (i in 1:number_of_sites) params[[paste0("MODEL", i)]] <- as.numeric(model_finer[i])
    
    DVRparams <- rbind(DVRparams, as.data.frame(params, stringsAsFactors = FALSE))
    
  } else {
    # if DTH contained NA, just skip (or you can store NA rows)
    cat("Skipping taxa (NA DTH):", as.character(TaxaName[[1]]), "\n")
  }
  
  toc()
}

# ============================================================
# Export results
# ============================================================
out_path <- paste0("results/Calibration/Coarse-Fine/", experimentID, ".xlsx")

write_xlsx(
  list(
    Coarse_Exact_Minima = all_exact_minima_coarse,
    Fine_Exact_Minima   = all_exact_minima_fine,
    Finer_Exact_Minima  = all_exact_minima_finer,
    Best_Params_Finer   = DVRparams
  ),
  out_path
)

cat("✅ Saved:", out_path, "\n")
