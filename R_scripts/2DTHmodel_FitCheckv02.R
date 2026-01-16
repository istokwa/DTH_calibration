# ============================================================
# 📊 MODEL FIT CHECK — Single combined plot + 4 metric lines
# ============================================================

library(ggplot2)
library(dplyr)
library(readxl)
library(Metrics)
library(grid)   # for unit()

# --- Inputs ---
excel_file <- "data/processed/(12-22-2025)-DVIstartest3.xlsx"
sheet_name <- "Finer_Best_Params"
experimentID <- "DVIstartest3"

obs_cols  <- c("DTH1", "DTH2", "DTH3")
pred_cols <- c("MODEL1", "MODEL2", "MODEL3")
sites     <- c("ISG1", "ISG2", "NGY1")

output_dir <- "results/1 - Regression/ModelFitCheck_Original (Test)"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# --- Read data ---
df <- read_excel(excel_file, sheet = sheet_name)

# --- Long format ---
combined_data <- bind_rows(lapply(1:3, function(i) {
  data.frame(
    Observed  = df[[obs_cols[i]]],
    Predicted = df[[pred_cols[i]]],
    Site      = sites[i]
  )
})) %>%
  filter(complete.cases(Observed, Predicted)) %>%
  mutate(Site = factor(Site, levels = sites))

# --- Metric function ---
calc_metrics <- function(obs, pred) {
  lm_fit <- lm(pred ~ obs)
  r2   <- summary(lm_fit)$r.squared
  rmse_v <- rmse(obs, pred)
  nse  <- 1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2)
  data.frame(R2 = r2, RMSE = rmse_v, NSE = nse)
}

# --- Per-site metrics ---
metrics_site <- combined_data %>%
  group_by(Site) %>%
  summarise(
    R2   = calc_metrics(Observed, Predicted)$R2,
    RMSE = calc_metrics(Observed, Predicted)$RMSE,
    NSE  = calc_metrics(Observed, Predicted)$NSE,
    .groups = "drop"
  )

# --- All-sites metrics ---
metrics_all <- calc_metrics(combined_data$Observed, combined_data$Predicted) %>%
  mutate(Site = "All sites")

# --- Build label table (4 lines) ---
metrics_df <- bind_rows(metrics_site, metrics_all) %>%
  mutate(
    Site = factor(Site, levels = c(sites, "All sites")),
    label = paste0(
      Site, ": R²=", sprintf("%.3f", R2),
      ", RMSE=", sprintf("%.3f", RMSE),
      ", NSE=", sprintf("%.3f", NSE)
    )
  )

# --- Place labels (top-left, stacked) ---
x_min <- min(combined_data$Observed, na.rm = TRUE)
x_max <- max(combined_data$Observed, na.rm = TRUE)
y_max <- max(combined_data$Predicted, na.rm = TRUE)

metrics_df <- metrics_df %>%
  arrange(Site) %>%
  mutate(
    x = x_min + 0.02 * (x_max - x_min),   # a bit right from left edge
    y = y_max,                             # anchor at top
    vjust = seq(0, by = 1.25, length.out = n()) # stack downward
  )

# --- Plot ---
p <- ggplot(combined_data, aes(x = Observed, y = Predicted, color = Site)) +
  geom_point(size = 2.5, alpha = 0.85) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 1) +
  geom_text(
    data = metrics_df,
    aes(x = x, y = y, label = label, vjust = vjust),
    inherit.aes = FALSE,
    hjust = 0,
    size = 6,
    color = "black"
  ) +
  labs(
    title = "Observed vs Predicted DTH",
    subtitle = paste0("ExperimentID: ", experimentID),
    x = "Observed DTH",
    y = "Predicted DTH"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(size = 18, face = "bold"),
    axis.text  = element_text(size = 14),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 14),
    legend.key.size = unit(1.4, "lines"),
    plot.margin = margin(t = 30, r = 15, b = 10, l = 20),
    plot.subtitle = element_text(size = 16, hjust = 0.5)
  )

# --- Save ---
ggsave(
  filename = file.path(output_dir, paste0("Observed_vs_Predicted_AllSites_", experimentID, ".png")),
  plot = p, width = 7, height = 6, dpi = 300
)

p
