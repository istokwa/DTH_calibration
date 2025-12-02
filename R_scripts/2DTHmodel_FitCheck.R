# ============================================================
# 📊 MODEL FIT CHECK — Multi-sheet version with SAVE OPTION
# ============================================================

# Load necessary libraries
library(ggplot2)
library(Metrics)
library(dplyr)
library(readxl)

# --- 1️⃣ File and sheets ---
excel_file <- "results/DVRparams.xlsx"
target_sheets <- c("Sheet1")

# --- 2️⃣ Output directory (set your preferred path here) ---
output_dir <- "results/ModelFitCheck_Original"   # <-- change this path if needed
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
    Observed = c(df[[obs_cols[1]]], df[[obs_cols[2]]], df[[obs_cols[3]]]),
    Predicted = c(df[[pred_cols[1]]], df[[pred_cols[2]]], df[[pred_cols[3]]]),
    Site = factor(rep(sites, each = nrow(df)))
  )
  
  # Combined plot
  combined_plot <- ggplot(combined_data, aes(x = Observed, y = Predicted, color = Site)) +
    geom_point(size = 2) +
    geom_smooth(method = "lm", se = FALSE) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
    labs(title = paste("Observed vs Predicted DTH - All Sites (", sheet_name, ")", sep = ""),
         x = "Observed DTH", y = "Predicted DTH") +
    theme_minimal() +
    theme(legend.position = "bottom") +
    geom_text(data = metrics_df,
              aes(x = min(combined_data$Observed),
                  y = max(combined_data$Predicted) - (0:2) * 2,
                  label = paste0(Site, ": R²=", round(R2,3),
                                 ", RMSE=", round(RMSE,2),
                                 ", NSE=", round(NSE,3))),
              color = "black", hjust = 0, vjust = 1, size = 4, inherit.aes = FALSE)
  
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