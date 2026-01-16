# ============================================================
# FINAL FULL VERSION (ENVIRONMENT + PHENOLOGY EMBEDDED)
# Among-environment tests per Metric (ANOVA/Tukey or Kruskal/Dunn)
# + CLD letters + facet_wrap in custom order
#
# INPUT: phenologydata.xlsx with columns: env, Temp, DL
#   - env: 12 sites
#   - each env has ~201 datapoints (time series)
#
# PHENOLOGY INCLUDED AS:
#   - Duration_days (n rows per env; phenology window proxy)
#   - Optional: early/mid/late phase summaries (commented block)
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(car)
  library(FSA)          # dunnTest
  library(multcompView) # multcompLetters
  library(ggplot2)
  library(scales)
  library(grid)         # unit()
})

# ============================================================
# USER SETTINGS
# ============================================================

setwd("D:/NU/Repository/GP_DTH-Rice/")

excel_file  <- "data/raw/phenologydata.xlsx"   # <-- adjust path if needed1
output_dir  <- "results/2 - Variance/EnvPheno"
actfilename <- "(01-09-2026)-Env_Phenology_Embedded"
experimentID <- "EnvPhenologyEmbedded"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

save_plot <- function(plot_obj, filename, width = 14, height = 7, dpi = 300) {
  filepath <- file.path(output_dir, paste0(filename, ".png"))
  ggsave(filename = filepath, plot = plot_obj, width = width, height = height, dpi = dpi)
  cat("✅ Saved plot:", filepath, "\n")
}

save_csv <- function(df, filename) {
  filepath <- file.path(output_dir, paste0(filename, ".csv"))
  write.csv(df, filepath, row.names = FALSE)
  cat("✅ Saved table:", filepath, "\n")
}

# ============================================================
# STEP 1: Load + QC
# ============================================================

EnvData <- read_excel(excel_file)

required_cols <- c("env", "Temp", "DL")
missing_cols <- setdiff(required_cols, colnames(EnvData))
if (length(missing_cols) > 0) {
  stop("Missing required columns in Excel: ", paste(missing_cols, collapse = ", "))
}

EnvData <- EnvData %>%
  mutate(
    env = as.factor(env),
    Temp = as.numeric(Temp),
    DL   = as.numeric(DL)
  )

cat("\n--- BASIC CHECKS ---\n")
cat("N rows:", nrow(EnvData), "\n")
cat("Environments:", paste(levels(EnvData$env), collapse = ", "), "\n")
cat("Rows per environment:\n")
print(table(EnvData$env))

cat("\n--- MISSING VALUE CHECK ---\n")
print(colSums(is.na(EnvData[, required_cols])))

# Add an explicit time index within each environment (phenology timeline)
# (This is useful if later you want phase-based summaries.)
EnvData <- EnvData %>%
  group_by(env) %>%
  mutate(DayIndex = row_number()) %>%
  ungroup()

# ============================================================
# STEP 2: Environment summaries (phenology embedded)
# ============================================================
# Phenology is "embedded" because Temp/DL are measured across a phenology timeline,
# and Duration_days is included as an explicit phenology proxy.

EnvSummary <- EnvData %>%
  group_by(env) %>%
  summarise(
    # Temperature summaries
    Temp_mean = mean(Temp, na.rm = TRUE),
    Temp_sd   = sd(Temp, na.rm = TRUE),
    Temp_min  = min(Temp, na.rm = TRUE),
    Temp_max  = max(Temp, na.rm = TRUE),
    
    # Daylength summaries
    DL_mean = mean(DL, na.rm = TRUE),
    DL_sd   = sd(DL, na.rm = TRUE),
    DL_min  = min(DL, na.rm = TRUE),
    DL_max  = max(DL, na.rm = TRUE),
    
    # Phenology proxy (duration of the time series per env)
    Duration_days = n(),
    
    .groups = "drop"
  )

cat("\n--- ENV SUMMARY (head) ---\n")
print(head(EnvSummary))

# ============================================================
# OPTIONAL (COMMENTED): add early/mid/late phase summaries
# If you want phenology explicitly split into phases (still environment-only).
# Uncomment this whole block if needed.
# ============================================================
# phase_breaks <- c(0, 0.33, 0.66, 1.0)
# EnvData <- EnvData %>%
#   group_by(env) %>%
#   mutate(
#     RelDay = DayIndex / max(DayIndex),
#     Phase = cut(RelDay,
#                 breaks = phase_breaks,
#                 labels = c("Early", "Mid", "Late"),
#                 include.lowest = TRUE)
#   ) %>%
#   ungroup()
#
# PhaseSummary <- EnvData %>%
#   group_by(env, Phase) %>%
#   summarise(
#     Temp_mean_phase = mean(Temp, na.rm = TRUE),
#     DL_mean_phase   = mean(DL,   na.rm = TRUE),
#     .groups = "drop"
#   ) %>%
#   pivot_wider(
#     names_from = Phase,
#     values_from = c(Temp_mean_phase, DL_mean_phase),
#     names_sep = "_"
#   )
#
# EnvSummary <- EnvSummary %>%
#   left_join(PhaseSummary, by = "env")

# ============================================================
# STEP 3: Long format + enforce facet order (CRITICAL)
# ============================================================

facet_levels <- c(
  "Temp_mean", "Temp_sd", "Temp_min", "Temp_max",
  "DL_mean",   "DL_sd",   "DL_min",   "DL_max",
  "Duration_days"
  # If you enabled phase summaries, add them here too, e.g.:
  # ,"Temp_mean_phase_Early","Temp_mean_phase_Mid","Temp_mean_phase_Late",
  # ,"DL_mean_phase_Early","DL_mean_phase_Mid","DL_mean_phase_Late"
)

long_env <- EnvSummary %>%
  pivot_longer(
    cols = all_of(facet_levels),
    names_to  = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    env    = as.factor(env),
    Metric = factor(Metric, levels = facet_levels)
  )

cat("\n✅ Long format complete.\n")

# ============================================================
# STEP 4: Among-environment tests per Metric + Posthoc + CLD
# ============================================================

param_test_results <- data.frame(
  Metric = character(),
  Test_Type = character(),
  Statistic = numeric(),
  p_value = numeric(),
  Interpretation = character(),
  stringsAsFactors = FALSE
)

posthoc_results <- data.frame(
  Metric = character(),
  Test_Type = character(),
  Comparison = character(),
  p_adj = numeric(),
  stringsAsFactors = FALSE
)

cld_results <- data.frame(
  Metric = character(),
  env = character(),
  CLD = character(),
  stringsAsFactors = FALSE
)

metrics <- levels(long_env$Metric)

for (m in metrics) {
  
  dfp <- long_env %>%
    filter(Metric == m) %>%
    filter(!is.na(Value), !is.na(env)) %>%
    mutate(env = droplevels(env))
  
  if (n_distinct(dfp$Value) <= 1) next
  if (n_distinct(dfp$env) < 2) next
  
  # Assumptions: residual normality + homogeneity
  normality_p <- tryCatch({
    fit <- lm(Value ~ env, data = dfp)
    shapiro.test(residuals(fit))$p.value
  }, error = function(e) NA_real_)
  
  var_homo_p <- tryCatch({
    car::leveneTest(Value ~ env, data = dfp)$`Pr(>F)`[1]
  }, error = function(e) NA_real_)
  
  parametric_ok <- !is.na(normality_p) && !is.na(var_homo_p) &&
    normality_p > 0.05 && var_homo_p > 0.05
  
  if (parametric_ok) {
    # ---------- ANOVA ----------
    aov_fit <- aov(Value ~ env, data = dfp)
    aov_sum <- summary(aov_fit)[[1]]
    stat <- aov_sum[["F value"]][1]
    pval <- aov_sum[["Pr(>F)"]][1]
    test_type <- "ANOVA"
    
    # Posthoc Tukey
    tk <- TukeyHSD(aov_fit, which = "env")$env
    tk_df <- as.data.frame(tk)
    tk_df$Comparison <- rownames(tk_df)
    
    posthoc_results <- rbind(
      posthoc_results,
      data.frame(
        Metric = as.character(m),
        Test_Type = "TukeyHSD",
        Comparison = tk_df$Comparison,
        p_adj = tk_df$`p adj`,
        stringsAsFactors = FALSE
      )
    )
    
    # CLD from Tukey adjusted p-values
    pvec <- tk_df$`p adj`
    names(pvec) <- tk_df$Comparison
    letters_obj <- multcompView::multcompLetters(pvec, threshold = 0.05)
    
    cld_results <- rbind(
      cld_results,
      data.frame(
        Metric = as.character(m),
        env = names(letters_obj$Letters),
        CLD = letters_obj$Letters,
        stringsAsFactors = FALSE
      )
    )
    
  } else {
    # ---------- Kruskal-Wallis ----------
    kw <- kruskal.test(Value ~ env, data = dfp)
    stat <- as.numeric(kw$statistic)
    pval <- as.numeric(kw$p.value)
    test_type <- "Kruskal-Wallis"
    
    # Posthoc Dunn (BH)
    dunn <- FSA::dunnTest(Value ~ env, data = dfp, method = "bh")$res
    
    posthoc_results <- rbind(
      posthoc_results,
      data.frame(
        Metric = as.character(m),
        Test_Type = "Dunn-BH",
        Comparison = dunn$Comparison,
        p_adj = dunn$P.adj,
        stringsAsFactors = FALSE
      )
    )
    
    # CLD from Dunn adjusted p-values
    pvec <- dunn$P.adj
    names(pvec) <- gsub(" - ", "-", dunn$Comparison)
    letters_obj <- multcompView::multcompLetters(pvec, threshold = 0.05)
    
    cld_results <- rbind(
      cld_results,
      data.frame(
        Metric = as.character(m),
        env = names(letters_obj$Letters),
        CLD = letters_obj$Letters,
        stringsAsFactors = FALSE
      )
    )
  }
  
  interpretation <- ifelse(pval < 0.05, "Significant", "Not significant")
  
  param_test_results <- rbind(
    param_test_results,
    data.frame(
      Metric = as.character(m),
      Test_Type = test_type,
      Statistic = stat,
      p_value = pval,
      Interpretation = interpretation,
      stringsAsFactors = FALSE
    )
  )
}

# Enforce order (CRITICAL if used in facet layers)
param_test_results$Metric <- factor(param_test_results$Metric, levels = facet_levels)
posthoc_results$Metric    <- factor(posthoc_results$Metric,    levels = facet_levels)
cld_results$Metric        <- factor(cld_results$Metric,        levels = facet_levels)

cat("\n--- TEST RESULTS ---\n")
print(param_test_results)

cat("\n--- POSTHOC (first 30 rows) ---\n")
print(head(posthoc_results, 30))

cat("\n--- CLD LETTERS ---\n")
print(cld_results)

# ============================================================
# STEP 5: CLD positions at top of each facet
# ============================================================

cld_positions <- long_env %>%
  group_by(Metric) %>%
  summarise(
    y_top = max(Value, na.rm = TRUE),
    y_range = diff(range(Value, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    y_cld = y_top + 0.10 * y_range
  ) %>%
  mutate(Metric = factor(Metric, levels = facet_levels))

cld_plot_df <- cld_results %>%
  left_join(cld_positions, by = "Metric") %>%
  mutate(Metric = factor(Metric, levels = facet_levels))

# ============================================================
# STEP 6: Plot (same visual language as your parameter plot)
# ============================================================

metric_labels <- c(
  Temp_mean = "Temp mean",
  Temp_sd   = "Temp SD",
  Temp_min  = "Temp min",
  Temp_max  = "Temp max",
  DL_mean   = "Daylength mean",
  DL_sd     = "Daylength SD",
  DL_min    = "Daylength min",
  DL_max    = "Daylength max",
  Duration_days = "Duration (days)"
)

p <- ggplot(long_env, aes(x = env, y = Value)) +
  geom_boxplot(
    outlier.colour = "darkgrey",
    outlier.size = 1,
    width = 0.6,
    box.linewidth = 1,
    whisker.linewidth = 1,
    median.linewidth = 1.25
  ) +
  facet_wrap(
    ~ Metric,
    scales = "free_y",
    strip.position = "top",
    labeller = labeller(Metric = metric_labels)
  ) +
  scale_y_continuous(
    limits = function(x) range(x, na.rm = TRUE),
    breaks = function(x) pretty(x, n = 5),
    labels = scales::label_number(accuracy = 0.01),
    expand = expansion(mult = c(0.05, 0.12))
  ) +
  theme_classic() +
  labs(x = "Environment (site)", y = "Value") +
  theme(
    panel.spacing = unit(0.5, "lines"),
    plot.margin = margin(5, 5, 5, 5),
    
    strip.text = element_text(size = 14, face = "bold"),
    
    axis.text.x = element_text(size = 11, hjust = 0.5, angle = 45),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title.y = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.ticks.x = element_line(linewidth = 1.2, colour = "black"),
    axis.ticks.y = element_line(linewidth = 1.2, colour = "black"),
    axis.ticks.length = unit(4, "pt"),
    
    axis.line = element_blank(),
    
    strip.background = element_rect(
      fill = "grey",
      colour = "black",
      linewidth = 1.5
    ),
    panel.border = element_rect(
      fill = NA,
      colour = "black",
      linewidth = 1.5
    ),
    
    strip.placement = "outside",
    strip.switch.pad.grid = unit(0, "pt"),
    strip.switch.pad.wrap = unit(0, "pt")
  )

p_final <- p +
  geom_text(
    data = cld_plot_df,
    aes(x = env, y = y_cld, label = CLD),
    inherit.aes = FALSE,
    size = 4.5,
    fontface = "bold",
    vjust = 0
  ) +
  labs(
    title = "Among-Environment Variation in Temperature & Daylength (Phenology Embedded)",
    subtitle = paste0("ExperimentID: ", experimentID, " | Each env has a phenology timeline (~201 days)")
  ) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 13, hjust = 0.5)
  )

print(p_final)

# Save plot + tables
save_plot(p_final, actfilename, width = 14, height = 7, dpi = 300)
save_csv(param_test_results, paste0(actfilename, "_TestResults"))
save_csv(posthoc_results,    paste0(actfilename, "_Posthoc"))
save_csv(cld_results,        paste0(actfilename, "_CLD"))

cat("\n✅ DONE. Outputs saved in: ", output_dir, "\n")
