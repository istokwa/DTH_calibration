# ============================================================
# FINAL FULL VERSION (DTH ONLY)
# Among-family tests per Site (DTH1/DTH2/DTH3)
# + CLD letters + facet_wrap (ISG1/ISG2/NGY1)
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(car)
  library(multcompView)
  library(FSA)        # dunnTest
  library(ggplot2)
  library(scales)
  library(readxl)
})

setwd("D:/NU/Repository/GP_DTH-Rice/")
output_dir <- "results/2 - Variance/"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

actfilename  <- "(12-15-2025)-Var-DTH_G_ISG2"
experimentID <- "G=ISG2"

save_plot <- function(plot_obj, filename, width = 10, height = 6, dpi = 300) {
  filepath <- file.path(output_dir, paste0(filename, ".png"))
  ggsave(filename = filepath, plot = plot_obj, width = width, height = height, dpi = dpi)
  cat("✅ Saved:", filepath, "\n")
}

# ============================================================
# STEP 1: Data QC + Family mapping
# ============================================================

OptParams <- read_excel("data/processed/(12-04-2025)-DVRparams-G_DTH2.xlsx")

family_map <- c(
  "WNAM_02" = "02",
  "WNAM_29" = "29",
  "WNAM_31" = "31",
  "WNAM_35" = "35",
  "WNAM_39" = "39",
  "WNAM_72" = "72",
  "WNAM_73" = "73"
)

numeric_cols <- c("DTH1", "DTH2", "DTH3", "MSE")

OptParams <- OptParams %>%
  mutate(
    Family = substr(Taxa, 1, 7),
    Family = ifelse(Family %in% names(family_map), family_map[Family], Family)
  )

cat("\n--- MISSING VALUE CHECK (DTH + MSE) ---\n")
na_summary <- colSums(is.na(OptParams[numeric_cols]))
print(na_summary)

OptParams_complete <- OptParams %>%
  filter(if_all(all_of(numeric_cols), ~ !is.na(.)))

mse_threshold <- quantile(OptParams_complete$MSE, 0.95, na.rm = TRUE)
OptParams_complete <- OptParams_complete %>%
  mutate(FitQuality = ifelse(MSE > mse_threshold, "BadFit", "GoodFit"))

cat("\nMSE 95th percentile threshold:", round(mse_threshold, 3), "\n")

# Keep only needed columns for this figure/analysis
DTH_df <- OptParams_complete %>%
  select(Taxa, Family, DTH1, DTH2, DTH3) %>%
  drop_na() %>%
  as.data.frame()

# ============================================================
# STEP 2: Long format + enforce facet order (DTH only)
# ============================================================

facet_levels <- c("DTH1", "DTH2", "DTH3")

long_DTH <- DTH_df %>%
  pivot_longer(
    cols = c(DTH1, DTH2, DTH3),
    names_to  = "Site",
    values_to = "DTH"
  ) %>%
  mutate(
    Family = as.factor(Family),
    Site = factor(Site, levels = facet_levels, labels = c("ISG1", "ISG2", "NGY1"))
  )

cat("\n✅ QC + long format complete (DTH only).\n")

# ============================================================
# STEP 3: Among-Family tests per Site + Posthoc + CLD
# ============================================================

site_test_results <- data.frame(
  Site = character(),
  Test_Type = character(),
  Statistic = numeric(),
  p_value = numeric(),
  Interpretation = character(),
  stringsAsFactors = FALSE
)

posthoc_results <- data.frame(
  Site = character(),
  Test_Type = character(),
  Comparison = character(),
  p_adj = numeric(),
  stringsAsFactors = FALSE
)

cld_results <- data.frame(
  Site = character(),
  Family = character(),
  CLD = character(),
  stringsAsFactors = FALSE
)

sites_levels <- levels(long_DTH$Site)

for (site in sites_levels) {

  dfp <- long_DTH %>%
    filter(Site == site) %>%
    filter(!is.na(DTH), !is.na(Family)) %>%
    mutate(Family = droplevels(Family))

  if (n_distinct(dfp$DTH) <= 1) next
  if (n_distinct(dfp$Family) < 2) next

  # Assumptions: residual normality + homogeneity
  normality_p <- tryCatch({
    fit <- lm(DTH ~ Family, data = dfp)
    shapiro.test(residuals(fit))$p.value
  }, error = function(e) NA_real_)

  var_homo_p <- tryCatch({
    car::leveneTest(DTH ~ Family, data = dfp)$`Pr(>F)`[1]
  }, error = function(e) NA_real_)

  parametric_ok <- !is.na(normality_p) && !is.na(var_homo_p) &&
    normality_p > 0.05 && var_homo_p > 0.05

  if (parametric_ok) {
    # ---------- ANOVA ----------
    aov_fit <- aov(DTH ~ Family, data = dfp)
    aov_sum <- summary(aov_fit)[[1]]
    stat <- aov_sum[["F value"]][1]
    pval <- aov_sum[["Pr(>F)"]][1]
    test_type <- "ANOVA"

    # Posthoc Tukey
    tk <- TukeyHSD(aov_fit, which = "Family")$Family
    tk_df <- as.data.frame(tk)
    tk_df$Comparison <- rownames(tk_df)

    posthoc_results <- rbind(
      posthoc_results,
      data.frame(
        Site = as.character(site),
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
        Site = as.character(site),
        Family = names(letters_obj$Letters),
        CLD = letters_obj$Letters,
        stringsAsFactors = FALSE
      )
    )

  } else {
    # ---------- Kruskal-Wallis ----------
    kw <- kruskal.test(DTH ~ Family, data = dfp)
    stat <- as.numeric(kw$statistic)
    pval <- as.numeric(kw$p.value)
    test_type <- "Kruskal-Wallis"

    # Posthoc Dunn (BH)
    dunn <- FSA::dunnTest(DTH ~ Family, data = dfp, method = "bh")$res

    posthoc_results <- rbind(
      posthoc_results,
      data.frame(
        Site = as.character(site),
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
        Site = as.character(site),
        Family = names(letters_obj$Letters),
        CLD = letters_obj$Letters,
        stringsAsFactors = FALSE
      )
    )
  }

  interpretation <- ifelse(pval < 0.05, "Significant", "Not significant")

  site_test_results <- rbind(
    site_test_results,
    data.frame(
      Site = as.character(site),
      Test_Type = test_type,
      Statistic = stat,
      p_value = pval,
      Interpretation = interpretation,
      stringsAsFactors = FALSE
    )
  )
}

cat("\n--- FAMILY COMPARISON TEST RESULTS (DTH) ---\n")
print(site_test_results)

cat("\n--- POST-HOC RESULTS (first 30 rows) ---\n")
print(head(posthoc_results, 30))

cat("\n--- CLD LETTERS (DTH) ---\n")
print(cld_results)

# ============================================================
# STEP 4: CLD positions at top of each facet (by Site)
# ============================================================

cld_positions <- long_DTH %>%
  group_by(Site) %>%
  summarise(
    y_top = max(DTH, na.rm = TRUE),
    y_range = diff(range(DTH, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    y_cld = y_top + 0.10 * y_range
  )

cld_plot_df <- cld_results %>%
  left_join(cld_positions, by = "Site")

# ============================================================
# STEP 5: Plot (DTH facets by Site; CLD at top)
# ============================================================

p <- ggplot(long_DTH, aes(x = Family, y = DTH)) +
  geom_boxplot(
    outlier.colour = "darkgrey",
    outlier.size = 1,
    width = 0.75,
    box.linewidth = 1,
    whisker.linewidth = 1,
    median.linewidth = 1.25
  ) +
  facet_wrap(~ Site, scales = "free_y", strip.position = "top", ncol = 1) +
  scale_y_continuous(
    limits = function(x) range(x, na.rm = TRUE),
    breaks = function(x) pretty(x, n = 5),
    labels = scales::label_number(accuracy = 0.1),
    expand = expansion(mult = c(0.05, 0.12))
  ) +
  theme_classic() +
  labs(x = "WNAM Family", y = "Days to Heading (DTH)") +
  theme(
    panel.spacing = unit(0.5, "lines"),
    plot.margin = margin(5, 5, 5, 5),

    strip.text = element_text(size = 14, face = "bold"),

    axis.text.x = element_text(size = 13, hjust = 0.5),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 15, face = "bold", hjust = 0.5),
    axis.title.y = element_text(size = 15, face = "bold", hjust = 0.5),

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
    aes(x = Family, y = y_cld, label = CLD),
    inherit.aes = FALSE,
    size = 5,
    fontface = "bold",
    vjust = 0
  ) +
  labs(
    title = "Among-Family Variation in DTH",
#    subtitle = paste0("ExperimentID: ", experimentID)
  ) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
#    plot.subtitle = element_text(size = 16, hjust = 0.5)
  )

print(p_final)

# Save
save_plot(p_final, actfilename, width = 5, height = 7, dpi = 300)
