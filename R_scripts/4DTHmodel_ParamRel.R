# ============================================================
# 📊 Pairwise Parameter Relationship Plots (Trade-off Analysis)
# ============================================================

library(readxl)
library(ggplot2)
library(dplyr)
library(viridis)

# --- 1️⃣ File and sheets ---
excel_file <- "results/GridSearch_Minima_Coarse_AllTaxa.xlsx"
target_sheets <- c("Coarse_Exact_Minima", "Coarse_Near_Minima", "Coarse_Best_Params")

# --- 2️⃣ Output directory ---
output_dir <- "results/Parameter_Relationship_Plots"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# --- 3️⃣ Function to plot bivariate parameter relationships ---
plot_param_relationship <- function(df, xvar, yvar, colorvar, sheet_name) {
  ggplot(df, aes_string(x = xvar, y = yvar, color = colorvar)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_smooth(method = "loess", color = "black", se = FALSE, span = 0.9) +
    scale_color_viridis_c(option = "plasma", direction = -1) +
    labs(
      title = paste0("Bivariate Parameter Relationship: ", xvar, " vs ", yvar),
      subtitle = paste0("Sheet: ", sheet_name, " — Color indicates ", colorvar, " (model fit error)"),
      x = paste0(xvar, " parameter"),
      y = paste0(yvar, " parameter"),
      color = paste0(colorvar, " (Model Error)")
    ) +
    theme_minimal(base_size = 13)
}

# --- 4️⃣ Loop through each sheet ---
for (sheet_name in target_sheets) {
  message("Processing sheet: ", sheet_name)
  df <- read_excel(excel_file, sheet = sheet_name)
  
  # Clean column names
  names(df) <- trimws(names(df))
  names(df) <- gsub(" ", "_", names(df))
  
  # Check necessary columns
  required_cols <- c("G", "A", "B", "Lc", "MSE")
  missing <- setdiff(required_cols, names(df))
  if (length(missing) > 0) {
    stop(paste("❌ Missing columns in", sheet_name, ":", paste(missing, collapse = ", ")))
  }
  
  # Output subfolder
  sheet_dir <- file.path(output_dir, sheet_name)
  dir.create(sheet_dir, recursive = TRUE, showWarnings = FALSE)
  
  # --- Plot 1: G vs A (parameter trade-off) ---
  p1 <- plot_param_relationship(df, "G", "A", "MSE", sheet_name)
  ggsave(file.path(sheet_dir, paste0("G_vs_A_", sheet_name, ".png")),
         plot = p1, width = 7, height = 5, dpi = 300)
  
  # --- Plot 2: B vs Lc (parameter compensation surface) ---
  p2 <- plot_param_relationship(df, "B", "Lc", "MSE", sheet_name)
  ggsave(file.path(sheet_dir, paste0("B_vs_Lc_", sheet_name, ".png")),
         plot = p2, width = 7, height = 5, dpi = 300)
  
  # --- Plot 3: Highlight best-fit region (top 1%) ---
  df <- df %>%
    mutate(FitGroup = ifelse(MSE <= min(MSE) * 1.01, "Top 1% (Lowest Error)", "Others"))
  
  p3 <- ggplot(df, aes(x = G, y = A, color = FitGroup)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_smooth(method = "lm", color = "grey40", se = FALSE) +
    scale_color_manual(values = c("Top 1% (Lowest Error)" = "red", "Others" = "lightblue")) +
    labs(
      title = paste("Parameter Trade-off: G vs A —", sheet_name),
      subtitle = "Highlighting lowest-error parameter combinations (Top 1%)",
      x = "G parameter",
      y = "A parameter",
      color = "Fit Category"
    ) +
    theme_minimal(base_size = 13)
  
  ggsave(file.path(sheet_dir, paste0("Tradeoff_G_vs_A_", sheet_name, ".png")),
         plot = p3, width = 7, height = 5, dpi = 300)
}

message("✅ All parameter relationship plots saved in: ", normalizePath(output_dir))

