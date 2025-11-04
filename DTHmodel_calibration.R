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

sourceCpp("mse_v02.cpp")

phenotypes <- read_excel("phenotypes_DVIDVR.xlsx")

# Generate DTH
dth <- data.frame(phenotypes[, 1])
colnames(dth) <- "taxa"
number_of_sites <- (ncol(phenotypes) - 1) / 3
for (i in 1:number_of_sites) {
  oldcolnames <- colnames(dth)
  dth[, i + 1] <- lapply(phenotypes[, i * 3 + 1] - phenotypes[, i * 3], as.numeric)
  colnames(dth) <- c(oldcolnames, paste("dth", i))
}

# Generate Daylenghts
daylegnths <- matrix(data = 0, nrow = 151, ncol = number_of_sites)
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
  day0 <- as.Date(phenotypes[[1, i * 3]])
  for (j in 1:151) {
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

# Generate Temperature
temperatures <- NULL

for (i in 1:number_of_sites) {
  temperaturedata <- read_excel(paste0(phenotypes[1, i * 3 - 1], ".xlsx"))
  day0 <- phenotypes[[1, i * 3]]
  row_of_day0 <- which(temperaturedata$Date == as.Date(day0))
  temp <- temperaturedata[row_of_day0:(row_of_day0 + 150), 2]
  temperatures <- if (is.null(temperatures))
    temp
  else
    cbind(temperatures, temp)
}

# DTH calculation using crop model
MSE <- function(dth,
                temperatures,
                daylegnths,
                G,
                Th,
                Lc,
                A,
                B,
                dvistar) {
  number_of_sites <- (ncol(phenotypes) - 1) / 3
  model <- vector()
  for (i in 1:number_of_sites) {
    model[i] <- calculateDTH(temperatures[, i], daylengths[, i], G, Th, Lc, A, B, DVIstar)
  }
  mse <- mean((model - dth)^2)
  return(mse)
}

# Coarse Space
bigspace <- list()
counter <- 1
for (G in seq(60, 10, by = -10)) {
  for (Th in seq(15, 30, by = 3)) {
    for (Lc in seq(10, 17, by = 1)) {
      for (A in seq(0, 2, by = 0.1)) {
        for (B in seq(0, 10, by = 0.5)) {
          # Store the current combination in the list
          bigspace[[counter]] <- c(G, Th, Lc, A, B)
          counter <- counter + 1
        }
      }
    }
  }
}
bigspace <- as.data.frame(do.call(rbind, bigspace))
colnames(bigspace) <- c("G", "Th", "Lc", "A", "B")

# Fine Space
template_of_smallspace <- list()
counter <- 1
for (G in seq(10, -10, by = -2)) {
  for (Th in seq(-3, 3, by = 1)) {
    for (Lc in seq(-1, 1, by = 0.2)) {
      for (A in seq(-0.1, 0.1, by = 0.02)) {
        for (B in seq(-0.5, 0.5, by = 0.2)) {
          template_of_smallspace[[counter]] <- c(G, Th, Lc, A, B)
          counter <- counter + 1
        }
      }
    }
  }
}
template_of_smallspace <- as.data.frame(do.call(rbind, template_of_smallspace))
colnames(template_of_smallspace) <- c("G", "Th", "Lc", "A", "B")

# Finer Space
template_of_smallerspace <- list()
counter <- 1
for (G in seq(2, -2, by = -0.5)) {
  for (Th in seq(-1, 1, by = 0.25)) {
    for (Lc in seq(-0.2, 0.2, by = 0.1)) {
      for (A in seq(-0.02, 0.02, by = 0.005)) {
        for (B in seq(-0.2, 0.2, by = 0.02)) {
          template_of_smallerspace[[counter]] <- c(G, Th, Lc, A, B)
          counter <- counter + 1
        }
      }
    }
  }
}
template_of_smallerspace <- as.data.frame(do.call(rbind, template_of_smallerspace))
colnames(template_of_smallerspace) <- c("G", "Th", "Lc", "A", "B")

# DF for results
DVRparams <- data.frame(matrix(ncol = 8 + 2 * number_of_sites, nrow = 0))
params <- vector("list",length = 8 + 2 * number_of_sites)
columnnames <- c("Taxa", "G", "Th", "Lc", "A", "B", "DVIStar", "MSE")
for (i in 1:number_of_sites) {
  columnnames <- c(columnnames, paste0("DTH", i))
}
for (i in 1:number_of_sites) {
  columnnames <- c(columnnames, paste0("MODEL", i))
}
colnames(DVRparams) <- columnnames
names(params) <- columnnames

# --- Coarse scan ---
all_exact_minima_coarse <- data.frame()
all_near_minima_coarse  <- data.frame()
all_best_params_coarse  <- data.frame()

# --- Fine scan ---
all_exact_minima_fine <- data.frame()
all_near_minima_fine  <- data.frame()
all_best_params_fine  <- data.frame()

# --- Finer scan ---
all_exact_minima_finer <- data.frame()
all_near_minima_finer  <- data.frame()
all_best_params_finer  <- data.frame()

# Scan taxa by taxa
for (taxa in 1:10) {
  tic()
  TaxaName <- "xxx"
  TaxaName <- phenotypes[taxa, 1]
  tempdth <- dth[taxa, 2:(number_of_sites + 1)]
  newMSE <- 99999.0
  if (!any(is.na(tempdth))) {
    getanMSE <- function(par) {
      G <- par[1]
      Th <- par[2]
      Lc <- par[3]
      A <- par[4]
      B <- par[5]
      tempMSE <- MSE_cpp(
        do.call(cbind, tempdth),
        do.call(cbind, temperatures),
        do.call(cbind, daylengths),
        G,
        Th,
        Lc,
        A,
        B,
        0
      ) # DVIstar=0
      return(tempMSE)
    }
    
    # coarse scan for bigspace grid
    results <- as.numeric(apply(bigspace, 1, getanMSE))
    min_index <- which.min(results)
    temp_best <- bigspace[min_index, ]
    min_mse <- results[min_index]
    results_df <- data.frame(bigspace, MSE = results)
    
    equal_minima <- results_df[results_df$MSE == min_mse, ]
    near_minima  <- results_df[results_df$MSE <= min_mse * 1.01, ]
    
    equal_minima$Taxa <- as.character(TaxaName[[1]])
    near_minima$Taxa  <- as.character(TaxaName[[1]])
    near_minima$IsBest <- near_minima$MSE == min_mse
    best_row <- data.frame(temp_best,
                           MSE = min_mse,
                           Taxa = as.character(TaxaName[[1]]),
                           stringsAsFactors = FALSE)
    
    # append
    all_exact_minima_coarse <- rbind(all_exact_minima_coarse, equal_minima)
    all_near_minima_coarse  <- rbind(all_near_minima_coarse, near_minima)
    all_best_params_coarse  <- rbind(all_best_params_coarse, best_row)
    
    cat("\n============================================\n")
    cat("Taxa:", TaxaName[[1]], "\n")
    cat("  ▶ Coarse best parameters:\n"); print(best_row)
    cat("\n  Exact minima:", nrow(equal_minima))
    cat("\n  Near minima (≤1%):", nrow(near_minima))
    cat("\n============================================\n\n")
    
    # fine scan
    print(paste("fine scan for", TaxaName))
    smallspace <- template_of_smallspace + as.vector(temp_best)
    results <- as.numeric(apply(smallspace, 1, getanMSE))
    min_index <- which.min(results)
    temp_best <- smallspace[min_index, ]
    min_mse <- results[min_index]
    results_df <- data.frame(smallspace, MSE = results)
    
    equal_minima <- results_df[results_df$MSE == min_mse, ]
    near_minima  <- results_df[results_df$MSE <= min_mse * 1.01, ]
    
    equal_minima$Taxa <- as.character(TaxaName[[1]])
    near_minima$Taxa  <- as.character(TaxaName[[1]])
    near_minima$IsBest <- near_minima$MSE == min_mse
    best_row <- data.frame(temp_best,
                           MSE = min_mse,
                           Taxa = as.character(TaxaName[[1]]),
                           stringsAsFactors = FALSE)
    
    # append
    all_exact_minima_fine <- rbind(all_exact_minima_fine, equal_minima)
    all_near_minima_fine  <- rbind(all_near_minima_fine, near_minima)
    all_best_params_fine  <- rbind(all_best_params_fine, best_row)
    
    cat("\n--------------------------------------------\n")
    cat("Taxa:", TaxaName[[1]], "\n")
    cat("  ▶ Fine best parameters:\n"); print(best_row)
    cat("\n  Exact minima:", nrow(equal_minima))
    cat("\n  Near minima (≤1%):", nrow(near_minima))
    cat("\n--------------------------------------------\n\n")
    
    
    
    
    # finer scan
    print(paste("finer scan for", TaxaName))
    smallerspace <- template_of_smallerspace + as.vector(temp_best)
    results <- as.numeric(apply(smallerspace, 1, getanMSE))
    min_index <- which.min(results)
    temp_best <- smallerspace[min_index, ]
    min_mse <- results[min_index]
    results_df <- data.frame(smallerspace, MSE = results)
    
    equal_minima <- results_df[results_df$MSE == min_mse, ]
    near_minima  <- results_df[results_df$MSE <= min_mse * 1.01, ]
    
    equal_minima$Taxa <- as.character(TaxaName[[1]])
    near_minima$Taxa  <- as.character(TaxaName[[1]])
    near_minima$IsBest <- near_minima$MSE == min_mse
    best_row <- data.frame(temp_best,
                           MSE = min_mse,
                           Taxa = as.character(TaxaName[[1]]),
                           stringsAsFactors = FALSE)
    
    # append
    all_exact_minima_finer <- rbind(all_exact_minima_finer, equal_minima)
    all_near_minima_finer  <- rbind(all_near_minima_finer, near_minima)
    all_best_params_finer  <- rbind(all_best_params_finer, best_row)
    
    cat("\n--------------------------------------------\n")
    cat("Taxa:", TaxaName[[1]], "\n")
    cat("  ▶ Finer best parameters:\n"); print(best_row)
    cat("\n  Exact minima:", nrow(equal_minima))
    cat("\n  Near minima (≤1%):", nrow(near_minima))
    cat("\n--------------------------------------------\n\n")
    
    
    
    
    tunedresult<-optim(par=temp_best,fn=getanMSE)
    params_o1<-c(1,tunedresult$par,NA,tunedresult$value,rep(NA,6))
    tunedresult<-optim(par=tunedresult$par,fn=getanMSE)
    params_o2<-c(2,tunedresult$par,NA,tunedresult$value,rep(NA,6))
    tunedresult<-optim(par=tunedresult$par,fn=getanMSE)
    params_o3<-c(3,tunedresult$par,NA,tunedresult$value,rep(NA,6))
    
    optG <- as.numeric(temp_best[1])
    optTh <- as.numeric(temp_best[2])
    optLc <- as.numeric(temp_best[3])
    optA <- as.numeric(temp_best[4])
    optB <- as.numeric(temp_best[5])
    optDVIstar <- 0
    newMSE <- as.numeric(getanMSE(unlist(temp_best)))
    
    model <- vector()
    for (i in 1:number_of_sites) {
      model[i] <- calculateDTH_optimized(temperatures[, i],
                                         daylengths[, i],
                                         optG,
                                         optTh,
                                         optLc,
                                         optA,
                                         optB,
                                         0)
    }
    
  } else {
    # if DTH contained NA
    optG <- NA
    optTh <- NA
    optLc <- NA
    optA <- NA
    optB <- NA
    optDVIstar <- NA
    tempdth <- rep(NA, number_of_sites)
    model <- rep(NA, number_of_sites)
    newMSE <- NA
  }
  
  params[1] <- TaxaName
  params[2] <- optG
  params[3] <- optTh
  params[4] <- optLc
  params[5] <- optA
  params[6] <- optB
  params[7] <- optDVIstar
  params[8] <- newMSE
  for (i in 1:length(tempdth)) {
    params[8 + i] <- tempdth[i]
  }
  for (i in 1:length(model)) {
    params[8 + i + length(tempdth)] <- model[i]
  }
  print(data.frame(params))
  DVRparams[taxa, ] <- params
  toc()
}

write_xlsx(DVRparams, "DVRparams_try_v03.xlsx")
