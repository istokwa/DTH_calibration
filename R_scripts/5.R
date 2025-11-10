### ============================================================
### STEP 0: Load Required Libraries
### ============================================================

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(car)
library(stats)
library(DescTools)
library(FSA)
library(multcomp)
library(rstatix)
library(agricolae)
library(ggpubr)
library(lme4)
library(FactoMineR)
library(factoextra)
library(cluster)
library(GGally)
library(Hmisc)
library(reshape2)
library(ggcorrplot)
library(multcompView)
library(tibble)
library(broom)

setwd("D:/NU/Repository/GP_DTH-Rice/")

### ============================================================
### GLOBAL SETTINGS FOR PLOT SAVING
### ============================================================

output_dir <- "results/Figures_PerDataset"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

save_plot <- function(plot_obj, filename, width = 8, height = 6, dpi = 300) {
  filepath <- file.path(output_dir, paste0(filename, ".png"))
  ggsave(filename = filepath, plot = plot_obj, width = width, height = height, dpi = dpi)
  cat("✅ Saved:", filepath, "\n")
}

### ============================================================
### STEP 1: Data QC (wide)
### ============================================================

OptParams <- read_excel("results/GridSearch_Minima_Coarse_AllTaxa.xlsx")

numeric_cols_wide <- c("G","Th","Lc","A","B","MSE",
                       "DTH1","DTH2","DTH3","MODEL1","MODEL2","MODEL3")

cat("\n--- MISSING VALUE CHECK (wide) ---\n")
na_summary <- colSums(is.na(OptParams[numeric_cols_wide]))
print(na_summary)

OptParams_complete_wide <- OptParams %>%
  filter(if_all(all_of(numeric_cols_wide), ~ !is.na(.)))

mse_threshold <- quantile(OptParams_complete_wide$MSE, 0.95, na.rm = TRUE)
OptParams_complete_wide <- OptParams_complete_wide %>%
  mutate(FitQuality = ifelse(MSE > mse_threshold, "BadFit", "GoodFit"))

cat("\nMSE 95th percentile threshold:", round(mse_threshold, 3), "\n")

cat("\n✅ QC complete.\n")

### ============================================================
### STEP 1b: Reshape Data to Dataset-based format
### ============================================================

# --- Long DTH reshape (safe version) ---
dth_long <- OptParams_complete_wide %>%
  select(Taxa, G, Th, Lc, A, B, MSE, DTH1, DTH2, DTH3) %>%
  pivot_longer(
    cols = starts_with("DTH"),
    names_to   = "Dataset_raw",
    values_to  = "DTH"
  ) %>%
  mutate(Dataset = case_when(
    Dataset_raw == "DTH1" ~ "Ishigaki1",
    Dataset_raw == "DTH2" ~ "Ishigaki2",
    Dataset_raw == "DTH3" ~ "Nagoya1",
    TRUE ~ Dataset_raw
  )) %>%
  select(-Dataset_raw)

# --- Long MODEL reshape (safe version) ---
model_long <- OptParams_complete_wide %>%
  select(Taxa, MODEL1, MODEL2, MODEL3) %>%
  pivot_longer(
    cols = starts_with("MODEL"),
    names_to   = "Dataset_raw",
    values_to  = "MODEL"
  ) %>%
  mutate(Dataset = case_when(
    Dataset_raw == "MODEL1" ~ "Ishigaki1",
    Dataset_raw == "MODEL2" ~ "Ishigaki2",
    Dataset_raw == "MODEL3" ~ "Nagoya1",
    TRUE ~ Dataset_raw
  )) %>%
  select(-Dataset_raw)

# Join both long tables by Taxa × Dataset
OptParams_long <- dth_long %>%
  left_join(model_long, by = c("Taxa", "Dataset")) %>%
  mutate(Dataset = factor(Dataset,
                          levels = c("Ishigaki1", "Ishigaki2", "Nagoya1")))

numeric_cols_long <- c("G","Th","Lc","A","B","MSE","DTH","MODEL")

cat("\n--- LONG-FORM DATA CREATED ---\n")
print(head(OptParams_long))

OptParams_complete <- OptParams_long %>%
  filter(if_all(all_of(numeric_cols_long), ~ !is.na(.)))

datasets <- unique(OptParams_complete$Dataset)

# Dataset-level summary
dataset_summary <- OptParams_complete %>%
  group_by(Dataset) %>%
  summarise(n = n(),
            across(all_of(numeric_cols_long),
                   ~ paste0(round(mean(.x, na.rm = TRUE), 2),
                            " ± ",
                            round(sd(.x, na.rm = TRUE), 2))),
            .groups = "drop")

overall_summary <- OptParams_complete %>%
  summarise(n = n(),
            across(all_of(numeric_cols_long),
                   ~ paste0(round(mean(.x, na.rm = TRUE), 2),
                            " ± ",
                            round(sd(.x, na.rm = TRUE), 2)))) %>%
  mutate(Dataset = "Overall_AllDatasets") %>%
  select(Dataset, everything())

summary_table <- bind_rows(dataset_summary, overall_summary)
cat("\n--- DATASET-LEVEL SUMMARY TABLE ---\n")
print(summary_table)

cat("\n✅ Dataset-based summary complete.\n")

### ============================================================
### STEP 2: Among-Dataset Tests (ANOVA / Kruskal-Wallis)
### ============================================================

param_test_results <- data.frame(
  Variable = character(),
  Test_Type = character(),
  Statistic = numeric(),
  p_value = numeric(),
  Interpretation = character(),
  stringsAsFactors = FALSE
)

for (var in numeric_cols_long) {
  y <- OptParams_complete[[var]]
  group <- OptParams_complete$Dataset
  if (length(unique(y)) <= 1) next
  
  normality_p <- tryCatch({
    pvals <- tapply(y, group,
                    function(x) if (length(x) > 2) shapiro.test(x)$p.value else NA)
    min(pvals, na.rm = TRUE)
  }, error = function(e) NA)
  
  var_homo_p <- tryCatch({
    leveneTest(y ~ as.factor(group))$`Pr(>F)`[1]
  }, error = function(e) NA)
  
  parametric_ok <- !is.na(normality_p) && !is.na(var_homo_p) &&
    normality_p > 0.05 && var_homo_p > 0.05
  
  if (parametric_ok) {
    test_res <- aov(y ~ as.factor(group))
    stat <- summary(test_res)[[1]][["F value"]][1]
    pval <- summary(test_res)[[1]][["Pr(>F)"]][1]
    test_type <- "ANOVA"
  } else {
    test_res <- kruskal.test(y ~ as.factor(group))
    stat <- test_res$statistic
    pval <- test_res$p.value
    test_type <- "Kruskal-Wallis"
  }
  
  interpretation <- ifelse(pval < 0.05, "Significant", "Not significant")
  
  param_test_results <- rbind(param_test_results,
                              data.frame(Variable = var,
                                         Test_Type = test_type,
                                         Statistic = stat,
                                         p_value = pval,
                                         Interpretation = interpretation))
}
cat("\n--- DATASET COMPARISON TEST RESULTS ---\n")
print(param_test_results)

### ============================================================
### STEP 3: Within-Dataset Tests (Variation Among Taxa)
### ============================================================

within_dataset_results <- data.frame(
  Dataset = character(),
  Variable = character(),
  Test_Type = character(),
  Statistic = numeric(),
  p_value = numeric(),
  Interpretation = character(),
  stringsAsFactors = FALSE
)

for (ds in datasets) {
  subset_data <- OptParams_complete %>% filter(Dataset == ds)
  
  for (var in numeric_cols_long) {
    y <- subset_data[[var]]
    taxa <- subset_data$Taxa
    
    if (length(unique(y)) <= 1) {
      within_dataset_results <- rbind(within_dataset_results,
                                      data.frame(Dataset = ds,
                                                 Variable = var,
                                                 Test_Type = NA,
                                                 Statistic = NA,
                                                 p_value = NA,
                                                 Interpretation = "All values identical; test skipped",
                                                 stringsAsFactors = FALSE))
      next
    }
    
    shapiro_p <- tryCatch({ shapiro.test(y)$p.value }, error = function(e) NA)
    levene_p <- tryCatch({ car::leveneTest(y ~ as.factor(taxa))$`Pr(>F)`[1] },
                         error = function(e) NA)
    
    parametric_ok <- !is.na(shapiro_p) && !is.na(levene_p) &&
      shapiro_p > 0.05 && levene_p > 0.05
    
    if (parametric_ok && length(unique(taxa)) > 1) {
      aov_res <- aov(y ~ as.factor(taxa))
      stat <- summary(aov_res)[[1]][["F value"]][1]
      pval <- summary(aov_res)[[1]][["Pr(>F)"]][1]
      test_type <- "ANOVA"
    } else if (length(unique(taxa)) > 1) {
      kw_res <- kruskal.test(y ~ as.factor(taxa))
      stat <- kw_res$statistic
      pval <- kw_res$p.value
      test_type <- "Kruskal-Wallis"
    } else {
      stat <- NA
      pval <- NA
      test_type <- NA
    }
    
    interpretation <- ifelse(!is.na(pval) & pval < 0.05, "Significant", "Not Significant")
    
    within_dataset_results <- rbind(within_dataset_results,
                                    data.frame(Dataset = ds,
                                               Variable = var,
                                               Test_Type = test_type,
                                               Statistic = round(stat, 4),
                                               p_value = signif(pval, 4),
                                               Interpretation = interpretation,
                                               stringsAsFactors = FALSE))
  }
}

cat("\n--- WITHIN-DATASET VARIATION (among taxa) ---\n")
print(within_dataset_results)

### ============================================================
### STEP 4: Post-hoc Pairwise Results (Among Datasets)
### ============================================================

posthoc_results <- data.frame(
  Variable = character(),
  Comparison = character(),
  p_value = numeric(),
  Adjusted_P = numeric(),
  Method = character(),
  Interpretation = character(),
  stringsAsFactors = FALSE
)

for (i in 1:nrow(param_test_results)) {
  var <- param_test_results$Variable[i]
  test_type <- param_test_results$Test_Type[i]
  pval <- param_test_results$p_value[i]
  if (is.na(pval) || pval >= 0.05) next
  
  y <- OptParams_complete[[var]]
  group <- as.factor(OptParams_complete$Dataset)
  
  if (test_type == "ANOVA") {
    aov_model <- aov(y ~ group)
    tukey_res <- TukeyHSD(aov_model)
    tmp <- as.data.frame(tukey_res$group)
    tmp$Comparison <- rownames(tmp)
    tmp <- tmp %>%
      mutate(Variable = var,
             p_value = `p adj`,
             Adjusted_P = `p adj`,
             Method = "Tukey HSD",
             Interpretation = ifelse(p_value < 0.05, "Significant", "Not significant")) %>%
      select(Variable, Comparison, p_value, Adjusted_P, Method, Interpretation)
    posthoc_results <- rbind(posthoc_results, tmp)
    
  } else if (test_type == "Kruskal-Wallis") {
    dunn_res <- dunnTest(y ~ group, method = "bonferroni")$res
    tmp <- dunn_res %>%
      mutate(Variable = var,
             p_value = P.unadj,
             Adjusted_P = P.adj,
             Method = "Dunn (Bonferroni)",
             Interpretation = ifelse(p_value < 0.05, "Significant", "Not significant")) %>%
      select(Variable, Comparison, p_value, Adjusted_P, Method, Interpretation)
    posthoc_results <- rbind(posthoc_results, tmp)
  }
}

cat("\n--- POST-HOC COMPARISONS (among datasets) ---\n")
print(posthoc_results)

significant_posthoc <- posthoc_results %>%
  filter(Adjusted_P < 0.05)
cat("\n--- SIGNIFICANT POST-HOC RESULTS ONLY ---\n")
print(significant_posthoc)

### ============================================================
### STEP 5: Boxplots WITH CLD by Dataset
### ============================================================

for (var in numeric_cols_long) {
  
  if (length(unique(OptParams_complete[[var]])) <= 1) next
  
  y <- OptParams_complete[[var]]
  group <- as.factor(OptParams_complete$Dataset)
  test_type <- param_test_results$Test_Type[param_test_results$Variable == var]
  cld_df <- NULL
  
  if (!is.na(test_type)) {
    if (test_type == "ANOVA") {
      aov_model <- aov(y ~ group)
      tukey_res <- glht(aov_model, linfct = mcp(group = "Tukey"))
      cld_letters <- cld(tukey_res)
      cld_df <- data.frame(Dataset = names(cld_letters$mcletters$Letters),
                           CLD = cld_letters$mcletters$Letters,
                           stringsAsFactors = FALSE)
    } else if (test_type == "Kruskal-Wallis") {
      dunn_res <- dunnTest(y ~ group, method = "bonferroni")$res
      ds <- levels(group)
      pmat <- matrix(1, nrow = length(ds), ncol = length(ds),
                     dimnames = list(ds, ds))
      for (i in 1:nrow(dunn_res)) {
        comps <- strsplit(dunn_res$Comparison[i], " - ")[[1]]
        pmat[comps[1], comps[2]] <- dunn_res$P.adj[i]
        pmat[comps[2], comps[1]] <- dunn_res$P.adj[i]
      }
      cld_letters <- multcompLetters(pmat, threshold = 0.05)
      cld_df <- data.frame(Dataset = names(cld_letters$Letters),
                           CLD = cld_letters$Letters,
                           stringsAsFactors = FALSE)
    }
  }
  
  max_y <- max(y, na.rm = TRUE)
  y_range <- max_y - min(y, na.rm = TRUE)
  
  p <- ggplot(OptParams_complete, aes(x = Dataset, y = .data[[var]], fill = Dataset)) +
    geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 1.5) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste("Variable:", var),
         x = "Dataset", y = var)
  
  if (!is.null(cld_df)) {
    cld_df <- cld_df %>% mutate(ypos = max_y + 0.1 * y_range)
    p <- p + geom_text(data = cld_df,
                       aes(x = Dataset, y = ypos, label = CLD),
                       inherit.aes = FALSE, vjust = 0)
  }
  
  print(p)
  save_plot(p, paste0("CLD_", var, "_byDataset"))
}
cat("\n✅ CLD-based dataset plots saved.\n")

### ============================================================
### STEP 6: Phenotypic Correlations (Overall + per Dataset)
### ============================================================

numeric_params <- c("G", "Th", "Lc", "A", "B")

cor_sig_fn <- function(data, mapping, ...) {
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  cor_res <- tryCatch(Hmisc::rcorr(x, y, type = "pearson"), error = function(e) NULL)
  if (!is.null(cor_res)) {
    r_val <- round(cor_res$r[1,2], 2)
    p_val <- cor_res$P[1,2]
    sig_star <- ifelse(p_val < 0.001, "***",
                       ifelse(p_val < 0.01, "**",
                              ifelse(p_val < 0.05, "*", "")))
    ggplot(data = data, mapping = mapping) +
      annotate("text", x = 0.5, y = 0.5,
               label = paste0(r_val, sig_star),
               size = 5,
               color = ifelse(r_val > 0.5, "red",
                              ifelse(r_val < -0.5, "blue", "black"))) +
      theme_void()
  } else { ggplot() + theme_void() }
}

subset_all <- OptParams_complete %>%
  select(all_of(numeric_params)) %>%
  drop_na()

if (nrow(subset_all) >= 2) {
  ggp_all <- ggpairs(subset_all,
                     lower = list(continuous = wrap("points", alpha = 0.6)),
                     diag  = list(continuous = "densityDiag"),
                     upper = list(continuous = cor_sig_fn)) +
    ggtitle("Phenotypic Correlations (Parameters) - All Datasets Pooled")
  print(ggp_all)
  save_plot(ggp_all, "Correlation_AllDatasets")
}

for (ds in datasets) {
  subset_ds <- OptParams_complete %>%
    filter(Dataset == ds) %>%
    select(all_of(numeric_params)) %>%
    drop_na()
  if (nrow(subset_ds) < 2) {
    cat("Skipping", ds, "- not enough taxa for correlation\n")
    next
  }
  ggp_ds <- ggpairs(subset_ds,
                    lower = list(continuous = wrap("points", alpha = 0.6)),
                    diag  = list(continuous = "densityDiag"),
                    upper = list(continuous = cor_sig_fn)) +
    ggtitle(paste("Phenotypic Correlations - Dataset:", ds))
  print(ggp_ds)
  save_plot(ggp_ds, paste0("Correlation_", ds))
}

### ============================================================
### STEP 7: Dataset-wise G×A and B×Lc Correlations (Table)
### ============================================================

compute_group_correlation <- function(data, var1, var2) {
  cor_test <- tryCatch(cor.test(data[[var1]], data[[var2]],
                                use = "complete.obs", method = "pearson"),
                       error = function(e) NULL)
  if (!is.null(cor_test)) {
    r_val <- cor_test$estimate
    p_val <- cor_test$p.value
  } else {
    r_val <- NA
    p_val <- NA
  }
  
  strength <- dplyr::case_when(
    abs(r_val) >= 0.9 ~ "Very strong",
    abs(r_val) >= 0.7 ~ "Strong",
    abs(r_val) >= 0.5 ~ "Moderate",
    abs(r_val) >= 0.3 ~ "Weak",
    TRUE ~ "Negligible"
  )
  
  direction <- ifelse(r_val > 0, "Positive", "Negative")
  
  significance <- dplyr::case_when(
    p_val < 0.001 ~ "***",
    p_val < 0.01  ~ "**",
    p_val < 0.05  ~ "*",
    TRUE ~ "ns"
  )
  
  interpretation <- paste(strength, direction, "correlation (", significance, ")", sep = " ")
  
  return(data.frame(r = r_val, p = p_val, Strength = strength,
                    Direction = direction, Significance = significance,
                    Interpretation = interpretation))
}

results_GA <- data.frame()
results_BLc <- data.frame()

for (ds in datasets) {
  ds_data <- OptParams_complete %>% filter(Dataset == ds)
  if (nrow(ds_data) > 3) {
    res_GA <- compute_group_correlation(ds_data, "G", "A") %>%
      mutate(Dataset = ds, n_Taxa = nrow(ds_data), Pair = "G × A")
    res_BLc <- compute_group_correlation(ds_data, "B", "Lc") %>%
      mutate(Dataset = ds, n_Taxa = nrow(ds_data), Pair = "B × Lc")
    results_GA <- bind_rows(results_GA, res_GA)
    results_BLc <- bind_rows(results_BLc, res_BLc)
  }
}

correlation_summary <- bind_rows(results_GA, results_BLc) %>%
  select(Dataset, Pair, n_Taxa, r, p, Strength, Direction, Significance, Interpretation) %>%
  arrange(Pair, desc(abs(r)))

cat("\n--- DATASET-WISE G×A and B×Lc CORRELATIONS ---\n")
print(correlation_summary)

### ============================================================
### STEP 8a: Dataset-level (Genetic Proxy) Correlation [Final Safe Version]
### ============================================================

param_cols <- c("G", "Th", "Lc", "A", "B")

dataset_means <- OptParams_complete %>%
  group_by(Dataset) %>%
  summarise(across(all_of(param_cols), \(x) mean(x, na.rm = TRUE)))

if (nrow(dataset_means) >= 2) {
  suppressWarnings({
    genetic_cor <- cor(dataset_means[, param_cols],
                       use = "pairwise.complete.obs",
                       method = "pearson")
  })
  
  # Check if all entries are NA
  if (all(is.na(genetic_cor))) {
    cat("\n⚠️  All correlations are NA — likely no between-dataset variation in parameters.\n")
  } else {
    cat("\n--- Dataset-level Parameter Correlations (Proxy for Genetic Component) ---\n")
    print(round(genetic_cor, 3))
    
    # Replace any remaining NAs with 0 to allow plotting
    genetic_cor[is.na(genetic_cor)] <- 0
    
    p_genetic <- tryCatch({
      ggcorrplot(
        genetic_cor,
        hc.order = TRUE,
        type = "lower",
        lab = TRUE,
        lab_size = 4,
        title = "Dataset-level Correlation (Genetic Proxy)",
        outline.color = "white",
        colors = c("blue", "white", "red")
      )
    }, error = function(e) {
      cat("⚠️  Could not generate correlation heatmap (", e$message, ")\n")
      NULL
    })
    
    if (!is.null(p_genetic)) {
      print(p_genetic)
      save_plot(p_genetic, "GeneticProxy_Correlation")
    }
  }
  
} else {
  cat("\n⚠️  Not enough datasets to compute correlation matrix.\n")
}

### ============================================================
### STEP 8b: Within-Dataset (Residual) Correlations  —  Safe
### ============================================================

within_dataset_cor_list <- list()

for (ds in unique(OptParams_complete$Dataset)) {
  subset_data <- OptParams_complete %>%
    filter(Dataset == ds) %>%
    select(all_of(param_cols)) %>%
    drop_na()
  
  if (nrow(subset_data) > 3) {
    suppressWarnings({
      cor_res <- cor(subset_data, use = "pairwise.complete.obs", method = "pearson")
    })
    within_dataset_cor_list[[ds]] <- cor_res
  }
}

if (length(within_dataset_cor_list) > 0) {
  within_dataset_summary <- do.call(rbind, lapply(names(within_dataset_cor_list), function(ds) {
    cor_mat <- within_dataset_cor_list[[ds]]
    data.frame(
      Dataset = ds,
      expand.grid(Var1 = rownames(cor_mat), Var2 = colnames(cor_mat)),
      Corr = as.vector(cor_mat)
    )
  })) %>%
    filter(Var1 != Var2)
  
  focus_pairs <- within_dataset_summary %>%
    filter((Var1 == "G" & Var2 == "A") | (Var1 == "B" & Var2 == "Lc"))
  
  cat("\n--- WITHIN-DATASET (RESIDUAL) CORRELATIONS: G×A and B×Lc ---\n")
  print(focus_pairs)
  
  avg_within_cor <- Reduce("+", within_dataset_cor_list) / length(within_dataset_cor_list)
  
  p_resid <- ggcorrplot(
    avg_within_cor,
    hc.order = TRUE,
    type = "lower",
    lab = TRUE,
    lab_size = 4,
    title = "Average Within-Dataset Correlation (Residual Component)",
    outline.color = "white",
    colors = c("blue", "white", "red")
  )
  print(p_resid)
  save_plot(p_resid, "Residual_Correlation")
  
} else {
  cat("\n⚠️  Not enough taxa per dataset to compute residual correlations.\n")
}


### ============================================================
### STEP 8c: Compare Correlation Types  —  Safe
### ============================================================

numeric_params <- c("G", "Th", "Lc", "A", "B")

# Phenotypic (pooled)
suppressWarnings({
  phenotypic_corr <- cor(OptParams_complete[, numeric_params], use = "pairwise.complete.obs")
})
r_GA_pheno  <- phenotypic_corr["G", "A"]
r_BLc_pheno <- phenotypic_corr["B", "Lc"]

# Dataset-level (proxy)
r_GA_genetic <- if (exists("genetic_cor")) genetic_cor["G", "A"] else NA
r_BLc_genetic <- if (exists("genetic_cor")) genetic_cor["B", "Lc"] else NA

# Residual (within-dataset mean)
if (length(within_dataset_cor_list) > 0) {
  within_corr_df <- do.call(rbind, lapply(within_dataset_cor_list, function(mat) {
    c(GA = mat["G", "A"], BLc = mat["B", "Lc"])
  }))
  r_GA_resid  <- mean(within_corr_df[, "GA"],  na.rm = TRUE)
  r_BLc_resid <- mean(within_corr_df[, "BLc"], na.rm = TRUE)
} else {
  r_GA_resid <- r_BLc_resid <- NA
}

corr_summary <- data.frame(
  Pair       = c("G×A", "B×Lc"),
  Phenotypic = c(r_GA_pheno,  r_BLc_pheno),
  Dataset    = c(r_GA_genetic, r_BLc_genetic),
  Residual   = c(r_GA_resid,  r_BLc_resid)
)

cat("\n--- SUMMARY: CORRELATION TYPE COMPARISON ---\n")
print(corr_summary %>% mutate(across(where(is.numeric), round, 3)))

corr_long <- corr_summary %>%
  pivot_longer(cols = c(Phenotypic, Dataset, Residual),
               names_to = "Correlation_Type",
               values_to = "r_value")

corr_type_plot <- ggplot(corr_long,
                         aes(x = Pair, y = r_value, fill = Correlation_Type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black") +
  geom_text(aes(label = ifelse(is.na(r_value), "NA", sprintf("%.2f", r_value))),
            position = position_dodge(width = 0.8),
            vjust = ifelse(corr_long$r_value > 0, -0.6, 1.3),
            size = 4.5) +
  scale_fill_manual(values = c("Phenotypic"="#6CA0DC",
                               "Dataset"="#E15759",
                               "Residual"="#59A14F")) +
  labs(title = "Comparison of Correlation Types (Phenotypic vs Dataset vs Residual)",
       x = "Parameter Pair", y = "Correlation Coefficient (r)",
       fill = "Correlation Type") +
  geom_hline(yintercept = 0, linetype="dashed", color="gray40") +
  theme_minimal(base_size = 14) +
  theme(legend.position="top",
        panel.grid.minor = element_blank(),
        plot.title = element_text(face="bold", size=16))
print(corr_type_plot)
save_plot(corr_type_plot, "CorrelationType_Comparison_DatasetLevel")


### ============================================================
### STEP 9: Relate to Performance (Observed DTH)
### ============================================================

dth_perf <- OptParams_complete %>%
  group_by(Dataset) %>%
  summarise(across(c(G, Th, Lc, A, B, DTH, MODEL), mean, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(DTH_diff = DTH - MODEL)

param_model <- lm(DTH ~ G + Th + Lc + A + B, data = dth_perf)

cat("\n--- Regression Summary: Dataset Mean DTH ~ Parameters ---\n")
print(summary(param_model))

r2_param <- summary(param_model)$r.squared
cat("\nR² (Proportion of variance explained by parameters):", round(r2_param, 3), "\n")

p_plot <- ggplot(dth_perf, aes(x = MODEL, y = DTH, color = Dataset)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(title = "Observed vs Predicted DTH per Dataset",
       x = "Predicted (MODEL)", y = "Observed (DTH)")
print(p_plot)
save_plot(p_plot, "DTH_vs_MODEL_byDataset")

cat("\n✅ Correlation and regression integration complete.\n")
