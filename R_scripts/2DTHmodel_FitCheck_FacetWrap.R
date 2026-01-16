library(ggplot2)
library(dplyr)
library(readxl)
library(purrr)
library(Metrics)
library(tibble)
library(ggrepel)

# ----------------------------
# 1) Define your experiments (edit paths)
# ----------------------------
experiments <- tibble(
  Experiment = c("BASE", "G=ISG2", "G=minDTH", "DVIstar"),
  File = c(
    "data/processed/(12-15-2025)-DVRparams-BASE.xlsx",
    "data/processed/(12-03-2025)-DVRparams-G_DTH2.xlsx",
    "data/processed/(12-08-2025)-DVRparams-G_lowestDTH.xlsx",
    "data/processed/(12-22-2025)-DVIstartest.xlsx"
  ),
  Sheet = c("Finer_Best_Params", "Finer_Best_Params", "Finer_Best_Params", "Finer_Best_Params")
)

# ----------------------------
# 2) Columns / sites
# ----------------------------
obs_cols  <- c("DTH1", "DTH2", "DTH3")
pred_cols <- c("MODEL1", "MODEL2", "MODEL3")
sites     <- c("ISG1", "ISG2", "NGY1")

output_dir <- "results/1 - Regression/FacetWrap_Experiments_Styled"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# 3) Read all experiments -> long data
# ----------------------------
read_one_experiment <- function(file, sheet, experiment_name) {
  df <- read_excel(file, sheet = sheet)
  
  bind_rows(lapply(1:3, function(i) {
    data.frame(
      Observed   = df[[obs_cols[i]]],
      Predicted  = df[[pred_cols[i]]],
      Site       = sites[i],
      Experiment = experiment_name
    )
  })) %>%
    filter(complete.cases(Observed, Predicted))
}

plot_data <- pmap_dfr(
  experiments,
  ~ read_one_experiment(file = ..2, sheet = ..3, experiment_name = ..1)
) %>%
  mutate(
    Site = factor(Site, levels = sites),
    Experiment = factor(Experiment, levels = experiments$Experiment)
  )

# ----------------------------
# 4) Metrics per experiment (site + all sites)
# ----------------------------
calc_metrics <- function(obs, pred) {
  lm_fit <- lm(pred ~ obs)
  r2   <- summary(lm_fit)$r.squared
  rmse <- Metrics::rmse(obs, pred)
  nse  <- 1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2)
  c(R2 = r2, RMSE = rmse, NSE = nse)
}

metrics_site <- plot_data %>%
  group_by(Experiment, Site) %>%
  summarise(
    R2   = calc_metrics(Observed, Predicted)["R2"],
    RMSE = calc_metrics(Observed, Predicted)["RMSE"],
    NSE  = calc_metrics(Observed, Predicted)["NSE"],
    .groups = "drop"
  ) %>%
  mutate(Site = as.character(Site))

metrics_all <- plot_data %>%
  group_by(Experiment) %>%
  summarise(
    R2   = calc_metrics(Observed, Predicted)["R2"],
    RMSE = calc_metrics(Observed, Predicted)["RMSE"],
    NSE  = calc_metrics(Observed, Predicted)["NSE"],
    .groups = "drop"
  ) %>%
  mutate(Site = "All sites")

metrics_df <- bind_rows(metrics_site, metrics_all) %>%
  mutate(
    line = paste0(
      Site, ": R²=", sprintf("%.3f", R2),
      ", RMSE=", sprintf("%.3f", RMSE),
      ", NSE=", sprintf("%.3f", NSE)
    )
  )

# One text block per Experiment (4 lines)
metrics_block <- metrics_df %>%
  group_by(Experiment) %>%
  summarise(label = paste(line, collapse = "\n"), .groups = "drop")

# ----------------------------
# 5) Put stats INSIDE each facet (top-left)
#    Use Inf/-Inf so it anchors consistently
# ----------------------------
metrics_block <- metrics_block %>%
  mutate(
    x = -Inf,   # left
    y = Inf     # top
  )

# ----------------------------
# 6) Plot
# ----------------------------
p <- ggplot(plot_data, aes(x = Observed, y = Predicted, color = Site)) +
  geom_point(size = 1, alpha = 0.85) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1.0) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 1.0) +
  facet_wrap(~ Experiment, ncol = 1) +
#  geom_text(
#    data = metrics_block,
#    aes(x = x, y = y, label = label),
#    inherit.aes = FALSE,
#    hjust = -0.02,  # nudge inside
#    vjust = 1.1,    # nudge down
#    size = 5.0,
#    color = "black",
#    lineheight = 1.1
#  ) +
  labs(
    title = "Observed vs Predicted DTH",
    #subtitle = "Facet by Experiment",
    x = "Observed DTH",
    y = "Predicted DTH",
    color = "Site"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    #plot.subtitle = element_text(size = 16, hjust = 0.5),
    
    panel.spacing = unit(0.5, "lines"),
    plot.margin = margin(5, 5, 5, 5),
    
    strip.background = element_rect(
      fill = "grey",
      colour = "black",
      linewidth = 1.5
    ),
    
    strip.text = element_text(size = 18, face = "bold"),
    strip.placement = "outside",
    strip.switch.pad.grid = unit(0, "pt"),
    strip.switch.pad.wrap = unit(0, "pt"),
    

    
    panel.border = element_rect(
      fill = NA,
      colour = "black",
      linewidth = 1.5
    ),
    
    axis.line = element_blank(),
    axis.ticks.y = element_line(linewidth = 1.2, colour = "black"),
    
    axis.title.x = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title.y = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.text  = element_text(size = 16),
    
    legend.position = "bottom",
    legend.title = element_text(size = 20, face = "bold"),
    legend.text  = element_text(size = 18),
    legend.key.size = unit(1.3, "lines"),
    
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
  ) +
  guides(color = guide_legend(override.aes = list(size = 5, linewidth = 1.5)))
  coord_cartesian(clip = "off")  # allow text nudges

# ----------------------------
# 7) Save (small & square-ish)
# ----------------------------
ggsave(
  filename = file.path(output_dir, "(12-22-2025)-DVIstartest.png"),
  plot = p,
  width = 4.5, height = 8, dpi = 300
)

p
