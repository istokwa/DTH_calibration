# ============================================================
# 📊 MODEL FIT CHECK — Multi-sheet version with SAVE OPTION
# ============================================================

# Load necessary libraries
library(ggplot2)
library(Metrics)
library(dplyr)
library(readxl)

# --- 1️⃣ File and sheets ---
excel_file <- "data/processed/MultiTaxa-G_DTH2-(12-04-2025).xlsx"
target_sheets <- c("Sheet1")

# --- 2️⃣ Output directory (set your preferred path here) ---
output_dir <- "results/ModelFitCheck-MultiTaxa-G_DTH2"   # <-- change this path if needed
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# --- 3️⃣ Column names and site labels ---
obs_cols <- c("DTH1", "DTH2", "DTH3")
pred_cols <- c("MODEL1", "MODEL2", "MODEL3")
sites <- c("Ishigaki1", "Ishigaki2", "Nagoya1")

# --- 4️⃣ Define plot function ---
plot_site <- function(obs, pred, site_name) {
  data <- na.omit(data.frame(Observed = obs, Predicted = pred))
  
  lm_fit <- lm(Predicted ~ Observed, data = data)
  r_squared <- summary(lm_fit)$r.squared
  rmse_val <- rmse(data$Observed, data$Predicted)
  nse_val <- 1 - sum((data$Observed - data$Predicted)^2) / sum((data$Observed - mean(data$Observed))^2)
  
  x_pos <- min(data$Observed) + 0.1 * diff(range(data$Observed))
  y_pos <- max(data$Predicted) - 0.1 * diff(range(data$Predicted))
  
  ggplot(data, aes(x = Observed, y = Predicted)) +
    geom_point(color = "blue", size = 3) +
    geom_smooth(method = "lm", color = "red", se = FALSE) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
    geom_text(aes(x = x_pos, y = y_pos,
                  label = paste0("R² = ", round(r_squared, 3),
                                 "\nRMSE = ", round(rmse_val, 2),
                                 "\nNSE = ", round(nse_val, 3))),
              hjust = 0, vjust = 1, color = "black", size = 5) +
    labs(title = paste("Observed vs Predicted DTH -", site_name),
         x = "Observed DTH", y = "Predicted DTH") +
    theme_minimal()
}

# --- 5️⃣ Function to process each sheet and save plots ---
analyze_sheet <- function(sheet_name) {
  message("Processing sheet: ", sheet_name)
  df <- read_excel(excel_file, sheet = sheet_name)
  
  # Output subfolder for this sheet
  sheet_dir <- file.path(output_dir, sheet_name)
  dir.create(sheet_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Per-site plots
  plot_list <- list()
  metrics_df <- data.frame(Site = sites, R2 = NA, RMSE = NA, NSE = NA)
  
  for (i in 1:3) {
    obs <- df[[obs_cols[i]]]
    pred <- df[[pred_cols[i]]]
    plot_obj <- plot_site(obs, pred, sites[i])
    
    # Save each site plot
    ggsave(filename = file.path(sheet_dir, paste0("Observed_vs_Predicted_", sites[i], ".png")),
           plot = plot_obj, width = 7, height = 5, dpi = 300)
    
    # Compute metrics
    valid <- complete.cases(obs, pred)
    obs <- obs[valid]; pred <- pred[valid]
    lm_fit <- lm(pred ~ obs)
    r_squared <- summary(lm_fit)$r.squared
    rmse_val <- rmse(obs, pred)
    nse_val <- 1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2)
    metrics_df[i, 2:4] <- c(r_squared, rmse_val, nse_val)
  }
  
  # Combined data for all sites
  combined_data <- data.frame(
    Observed  = c(df[[obs_cols[1]]], df[[obs_cols[2]]], df[[obs_cols[3]]]),
    Predicted = c(df[[pred_cols[1]]], df[[pred_cols[2]]], df[[pred_cols[3]]]),
    Site      = factor(rep(sites, each = nrow(df)))
  )
  
  ## --- Fit regressions for each site (for equations) ---
  site_levels <- levels(combined_data$Site)
  eq_df <- data.frame(
    Site      = site_levels,
    intercept = NA_real_,
    slope     = NA_real_,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(site_levels)) {
    subdat <- subset(combined_data, Site == site_levels[i])
    subdat <- na.omit(subdat)
    if (nrow(subdat) > 1) {
      lm_i <- lm(Predicted ~ Observed, data = subdat)
      eq_df$intercept[i] <- coef(lm_i)[1]
      eq_df$slope[i]     <- coef(lm_i)[2]
    }
  }
  # --- Fit regressions for each site (for equations + metrics) ---
  site_levels <- levels(combined_data$Site)
  
  eq_df <- data.frame(
    Site      = site_levels,
    intercept = NA_real_,
    slope     = NA_real_,
    R2        = NA_real_,
    RMSE      = NA_real_,
    NSE       = NA_real_,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(site_levels)) {
    subdat <- subset(combined_data, Site == site_levels[i])
    subdat <- na.omit(subdat)
    
    if (nrow(subdat) > 1) {
      lm_i <- lm(Predicted ~ Observed, data = subdat)
      
      eq_df$intercept[i] <- coef(lm_i)[1]
      eq_df$slope[i]     <- coef(lm_i)[2]
      eq_df$R2[i]        <- summary(lm_i)$r.squared
      eq_df$RMSE[i]      <- rmse(subdat$Observed, subdat$Predicted)
      eq_df$NSE[i]       <- 1 - sum((subdat$Observed - subdat$Predicted)^2) /
        sum((subdat$Observed - mean(subdat$Observed))^2)
    }
  }
  
  # --- GLOBAL regression (ALL sites combined) ---
  combined_clean <- na.omit(combined_data)
  lm_all <- lm(Predicted ~ Observed, data = combined_clean)
  coef_all <- coef(lm_all)
  
  R2_all  <- summary(lm_all)$r.squared
  RMSE_all <- rmse(combined_clean$Observed, combined_clean$Predicted)
  NSE_all  <- 1 - sum((combined_clean$Observed - combined_clean$Predicted)^2) /
    sum((combined_clean$Observed - mean(combined_clean$Observed))^2)
  
  eq_all <- data.frame(
    Site      = "AllSites",
    intercept = coef_all[1],
    slope     = coef_all[2],
    R2        = R2_all,
    RMSE      = RMSE_all,
    NSE       = NSE_all,
    stringsAsFactors = FALSE
  )
  
  # Combine site-wise + global
  eq_labels <- rbind(eq_df, eq_all)
  
  # Equation labels
  eq_labels$Equation <- sprintf(
    "%s: y = %.2f + %.2f x",
    eq_labels$Site, eq_labels$intercept, eq_labels$slope
  )
  
  # Metric labels
  eq_labels$Metrics <- sprintf(
    "%s: R²=%.3f, RMSE=%.2f, NSE=%.3f",
    eq_labels$Site, eq_labels$R2, eq_labels$RMSE, eq_labels$NSE
  )
  
  # Positioning for text
  x_min <- min(combined_clean$Observed, na.rm = TRUE)
  y_max <- max(combined_clean$Predicted, na.rm = TRUE)
  y_range <- diff(range(combined_clean$Predicted, na.rm = TRUE))
  
  spacing <- 0.06 * y_range
  eq_labels$x_eq <- x_min
  eq_labels$y_eq <- y_max - (seq_len(nrow(eq_labels)) - 1) * spacing
  
  # Metrics on the right side
  x_max <- max(combined_clean$Observed, na.rm = TRUE)
  eq_labels$x_met <- x_max
  eq_labels$y_met <- y_max - (seq_len(nrow(eq_labels)) - 1) * spacing
  
  
  # Combined plot
  combined_plot <- ggplot(combined_data, aes(x = Observed, y = Predicted, color = Site)) +
    geom_point(size = 2) +
    
    # per-site LM lines
    geom_smooth(method = "lm", se = FALSE) +
    
    # global LM line
    geom_abline(slope = coef_all[2], intercept = coef_all[1],
                linetype = "dotted", linewidth = 0.9, color = "black") +
    
    # 1:1 line
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
    
    labs(title = paste("Observed vs Predicted DTH - All Sites (", sheet_name, ")", sep = ""),
         x = "Observed DTH", y = "Predicted DTH") +
    theme_minimal() +
    theme(legend.position = "bottom") +
    
    # --- Regression equations (LEFT side) ---
    geom_text(data = eq_labels,
              aes(x = x_eq, y = y_eq, label = Equation),
              color = "black", hjust = 0, vjust = 1, size = 3.5, inherit.aes = FALSE) +
    
    # --- Metrics text (RIGHT side) ---
    geom_text(data = eq_labels,
              aes(x = x_met, y = y_met, label = Metrics),
              color = "black", hjust = 1, vjust = 1, size = 3.5, inherit.aes = FALSE)
  
  # Save combined plot
  ggsave(filename = file.path(sheet_dir, paste0("Combined_AllSites_", sheet_name, ".png")),
         plot = combined_plot, width = 8, height = 6, dpi = 300)
  
  # Save metrics table to CSV
  write.csv(metrics_df, file.path(sheet_dir, paste0("Metrics_", sheet_name, ".csv")), row.names = FALSE)
  
  # Return list for reference
  return(list(sheet = sheet_name, metrics = metrics_df, combined = combined_plot))
}

# --- 6️⃣ Run analysis for all sheets ---
results_list <- lapply(target_sheets, analyze_sheet)
names(results_list) <- target_sheets

message("✅ All plots and metrics saved in: ", normalizePath(output_dir))