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

output_dir <- "results/Figures_MultiBestMSE_G70"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

save_plot <- function(plot_obj, filename, width = 8, height = 6, dpi = 300) {
  filepath <- file.path(output_dir, paste0(filename, ".png"))
  ggsave(filename = filepath, plot = plot_obj, width = width, height = height, dpi = dpi)
  cat("✅ Saved:", filepath, "\n")
}

### ============================================================
### STEP 1: Data QC + Family-based Summary
### ============================================================

OptParams <- read_excel("results/Coarse_ExactMinima_MultiTaxaG70.xlsx")

family_map <- c(
  "WNAM_02" = "WNAM02",
  "WNAM_29" = "WNAM29",
  "WNAM_31" = "WNAM31",
  "WNAM_35" = "WNAM35",
  "WNAM_39" = "WNAM39",
  "WNAM_72" = "WNAM72",
  "WNAM_73" = "WNAM73"
)

OptParams <- OptParams %>%
  mutate(Family = substr(Taxa, 1, 7),
         Family = ifelse(Family %in% names(family_map),
                         family_map[Family],
                         Taxa))

numeric_cols <- c("G","Th","Lc","A","B","MSE",
                  "DTH1","DTH2","DTH3","MODEL1","MODEL2","MODEL3")

cat("\n--- MISSING VALUE CHECK ---\n")
na_summary <- colSums(is.na(OptParams[numeric_cols]))
print(na_summary)

OptParams_complete <- OptParams %>%
  filter(if_all(all_of(numeric_cols), ~ !is.na(.)))

mse_threshold <- quantile(OptParams_complete$MSE, 0.95, na.rm = TRUE)
OptParams_complete <- OptParams_complete %>%
  mutate(FitQuality = ifelse(MSE > mse_threshold, "BadFit", "GoodFit"))

cat("\nMSE 95th percentile threshold:", round(mse_threshold, 3), "\n")

family_summary <- OptParams_complete %>%
  group_by(Family) %>%
  summarise(n = n(),
            across(all_of(numeric_cols),
                   ~ paste0(round(mean(.x, na.rm = TRUE), 2),
                            " ± ",
                            round(sd(.x, na.rm = TRUE), 2))),
            .groups = "drop")

overall_summary <- OptParams_complete %>%
  summarise(n = n(),
            across(all_of(numeric_cols),
                   ~ paste0(round(mean(.x, na.rm = TRUE), 2),
                            " ± ",
                            round(sd(.x, na.rm = TRUE), 2)))) %>%
  mutate(Family = "Overall_NoGroup") %>%
  select(Family, everything())

summary_table <- bind_rows(family_summary, overall_summary)
cat("\n--- FAMILY-LEVEL SUMMARY TABLE ---\n")
print(summary_table)

cat("\n✅ QC and summary complete.\n")

### ============================================================
### STEP 2: Among-Family Tests (ANOVA / Kruskal-Wallis)
### ============================================================

param_test_results <- data.frame(
  Parameter = character(),
  Test_Type = character(),
  Statistic = numeric(),
  p_value = numeric(),
  Interpretation = character(),
  stringsAsFactors = FALSE
)

families <- unique(OptParams_complete$Family)

for (param in numeric_cols) {
  y <- OptParams_complete[[param]]
  group <- OptParams_complete$Family
  if (length(unique(y)) <= 1) next
  
  normality_p <- tryCatch({
    pvals <- tapply(y, group, function(x) if(length(x) > 2) shapiro.test(x)$p.value else NA)
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
                              data.frame(Parameter = param,
                                         Test_Type = test_type,
                                         Statistic = stat,
                                         p_value = pval,
                                         Interpretation = interpretation))
}
cat("\n--- FAMILY COMPARISON TEST RESULTS ---\n")
print(param_test_results)

### ============================================================
### STEP 3: Within-Family Tests (Optional, no plots)
### ============================================================

within_family_results <- data.frame(
  Family = character(),
  Parameter = character(),
  Test_Type = character(),
  Statistic = numeric(),
  p_value = numeric(),
  Interpretation = character(),
  stringsAsFactors = FALSE
)

for (fam in families) {
  subset_data <- OptParams_complete %>% filter(Family == fam)
  
  for (param in numeric_cols) {
    y <- subset_data[[param]]
    taxa <- subset_data$Taxa
    
    if (length(unique(y)) <= 1) {
      within_family_results <- rbind(within_family_results,
                                     data.frame(Family = fam,
                                                Parameter = param,
                                                Test_Type = NA,
                                                Statistic = NA,
                                                p_value = NA,
                                                Interpretation = "All values identical; test skipped",
                                                stringsAsFactors = FALSE))
      next
    }
    
    shapiro_p <- tryCatch({ shapiro.test(y)$p.value }, error = function(e) NA)
    levene_p <- tryCatch({ car::leveneTest(y ~ as.factor(taxa))$`Pr(>F)`[1] }, error = function(e) NA)
    
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
    
    within_family_results <- rbind(within_family_results,
                                   data.frame(Family = fam,
                                              Parameter = param,
                                              Test_Type = test_type,
                                              Statistic = round(stat, 4),
                                              p_value = signif(pval,4),
                                              Interpretation = interpretation,
                                              stringsAsFactors = FALSE))
  }
}

### ============================================================
### STEP 4: Post-hoc Pairwise Results (Tables only)
### ============================================================

posthoc_results <- data.frame(
  Parameter = character(),
  Comparison = character(),
  p_value = numeric(),
  Adjusted_P = numeric(),
  Method = character(),
  Interpretation = character(),
  stringsAsFactors = FALSE
)

for (i in 1:nrow(param_test_results)) {
  param <- param_test_results$Parameter[i]
  test_type <- param_test_results$Test_Type[i]
  pval <- param_test_results$p_value[i]
  if (is.na(pval) || pval >= 0.05) next
  
  y <- OptParams_complete[[param]]
  group <- as.factor(OptParams_complete$Family)
  
  if (test_type == "ANOVA") {
    aov_model <- aov(y ~ group)
    tukey_res <- TukeyHSD(aov_model)
    tmp <- as.data.frame(tukey_res$group)
    tmp$Comparison <- rownames(tmp)
    tmp <- tmp %>%
      mutate(Parameter = param,
             p_value = `p adj`,
             Adjusted_P = `p adj`,
             Method = "Tukey HSD",
             Interpretation = ifelse(p_value < 0.05, "Significant", "Not significant")) %>%
      select(Parameter, Comparison, p_value, Adjusted_P, Method, Interpretation)
    posthoc_results <- rbind(posthoc_results, tmp)
    
  } else if (test_type == "Kruskal-Wallis") {
    dunn_res <- dunnTest(y ~ group, method = "bonferroni")$res
    tmp <- dunn_res %>%
      mutate(Parameter = param,
             p_value = P.unadj,
             Adjusted_P = P.adj,
             Method = "Dunn (Bonferroni)",
             Interpretation = ifelse(p_value < 0.05, "Significant", "Not significant")) %>%
      select(Parameter, Comparison, p_value, Adjusted_P, Method, Interpretation)
    posthoc_results <- rbind(posthoc_results, tmp)
  }
}

cat("\n--- POST-HOC COMPARISONS ---\n")
print(posthoc_results)

significant_posthoc <- posthoc_results %>%
  filter(Adjusted_P < 0.05)
cat("\n--- SIGNIFICANT POST-HOC RESULTS ONLY ---\n")
print(significant_posthoc)

### ============================================================
### STEP 5: Boxplots WITH CLD (Compact Letter Display)
### ============================================================

for (param in numeric_cols) {
  
  if (length(unique(OptParams_complete[[param]])) <= 1) next
  
  y <- OptParams_complete[[param]]
  group <- as.factor(OptParams_complete$Family)
  test_type <- param_test_results$Test_Type[param_test_results$Parameter == param]
  cld_df <- NULL
  
  if (!is.na(test_type)) {
    if (test_type == "ANOVA") {
      aov_model <- aov(y ~ group)
      tukey_res <- glht(aov_model, linfct = mcp(group = "Tukey"))
      cld_letters <- cld(tukey_res)
      cld_df <- data.frame(Family = names(cld_letters$mcletters$Letters),
                           CLD = cld_letters$mcletters$Letters,
                           stringsAsFactors = FALSE)
    } else if (test_type == "Kruskal-Wallis") {
      dunn_res <- dunnTest(y ~ group, method = "bonferroni")$res
      fams <- unique(group)
      pmat <- matrix(1, nrow = length(fams), ncol = length(fams),
                     dimnames = list(fams, fams))
      for (i in 1:nrow(dunn_res)) {
        comps <- strsplit(dunn_res$Comparison[i], " - ")[[1]]
        pmat[comps[1], comps[2]] <- dunn_res$P.adj[i]
        pmat[comps[2], comps[1]] <- dunn_res$P.adj[i]
      }
      cld_letters <- multcompLetters(pmat, threshold = 0.05)
      cld_df <- data.frame(Family = names(cld_letters$Letters),
                           CLD = cld_letters$Letters,
                           stringsAsFactors = FALSE)
    }
  }
  
  max_y <- max(y, na.rm = TRUE)
  y_range <- max_y - min(y, na.rm = TRUE)
  
  p <- ggplot(OptParams_complete, aes(x = Family, y = .data[[param]], fill = Family)) +
    geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 1.5) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste("Parameter:", param),
         x = "Family", y = param)
  
  if (!is.null(cld_df)) {
    cld_df <- cld_df %>% mutate(ypos = max_y + 0.1 * y_range)
    p <- p + geom_text(data = cld_df,
                       aes(x = Family, y = ypos, label = CLD),
                       inherit.aes = FALSE, vjust = 0)
  }
  
  print(p)
  save_plot(p, paste0("CLD_", param))
}
cat("\n✅ CLD-based diagnostic plots saved.\n")

### ============================================================
### STEP 6: PCA + CLUSTERING
### ============================================================

param_matrix <- OptParams_complete %>%
  select(Taxa, G, Th, Lc, A, B) %>%
  drop_na() %>%
  as.data.frame()

# handle duplicate taxa
param_matrix$Taxa <- make.unique(as.character(param_matrix$Taxa))

pca_data <- OptParams_complete %>%
  select(Taxa, Family, G, Th, Lc, A, B) %>%
  drop_na() %>%
  mutate(Taxa = make.unique(as.character(Taxa)))

rownames(param_matrix) <- param_matrix$Taxa
param_matrix <- param_matrix[, -1]

pca_res <- FactoMineR::PCA(param_matrix, scale.unit = TRUE, graph = FALSE)

p_pca_ind <- fviz_pca_ind(
  pca_res,
  geom.ind = "point",
  col.ind = pca_data$Family,
  addEllipses = TRUE,
  palette = "Dark2",
  legend.title = "Family",
  repel = TRUE,
  title = "PCA of Optimized Parameters (G, Th, Lc, A, B)"
)
print(p_pca_ind)
save_plot(p_pca_ind, "PCA_Individuals")

p_pca_var <- fviz_pca_var(
  pca_res,
  col.var = "contrib",
  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
  title = "Contribution of Parameters to Principal Components"
)
print(p_pca_var)
save_plot(p_pca_var, "PCA_Variables")

family_means_pca <- pca_data %>%
  group_by(Family) %>%
  summarise(across(c(G, Th, Lc, A, B), mean, na.rm = TRUE)) %>%
  as.data.frame()
rownames(family_means_pca) <- family_means_pca$Family
family_means_pca <- family_means_pca[, -1]

dist_matrix <- dist(scale(family_means_pca))
hc <- hclust(dist_matrix, method = "ward.D2")
png(file.path(output_dir, "Hierarchical_Clustering.png"),
    width = 8, height = 6, units = "in", res = 300)
plot(hc, main = "Hierarchical Clustering of Families (Ward.D2)",
     xlab = "Family", sub = "", cex = 1.2)
dev.off()
cat("✅ Saved:", file.path(output_dir, "Hierarchical_Clustering.png"), "\n")

### ============================================================
### STEP 6b: PCA Loadings + Biological Interpretation
### ============================================================

pca_loadings <- as.data.frame(pca_res$var$coord)
pca_loadings$Parameter <- rownames(pca_loadings)

cat("\n--- PCA LOADINGS TABLE (rounded) ---\n")
print(round(pca_loadings[, sapply(pca_loadings, is.numeric)], 3))
cat("\nParameters:\n")
print(pca_loadings$Parameter)

loading_strength <- pca_loadings %>%
  mutate(across(starts_with("Dim."), abs)) %>%
  arrange(desc(Dim.1)) %>%
  select(Parameter, Dim.1, Dim.2, Dim.3, Dim.4, Dim.5)

cat("\n--- PARAMETERS CONTRIBUTING STRONGLY TO EACH DIMENSION ---\n")
loading_strength_rounded <- loading_strength
loading_strength_rounded[, sapply(loading_strength_rounded, is.numeric)] <-
  round(loading_strength_rounded[, sapply(loading_strength_rounded, is.numeric)], 3)
print(loading_strength_rounded)

top_contributors <- loading_strength %>%
  summarise(
    Dim1_Top = paste(head(Parameter[order(-Dim.1)], 3), collapse = ", "),
    Dim2_Top = paste(head(Parameter[order(-Dim.2)], 3), collapse = ", "),
    Dim3_Top = paste(head(Parameter[order(-Dim.3)], 3), collapse = ", "),
    Dim4_Top = paste(head(Parameter[order(-Dim.4)], 3), collapse = ", "),
    Dim5_Top = paste(head(Parameter[order(-Dim.5)], 3), collapse = ", ")
  )

cat("\n--- PCA INTERPRETATION SUMMARY ---\n")
cat(paste0(
  "Dim.1 is mainly driven by: ", top_contributors$Dim1_Top, "\n",
  "Dim.2 is mainly driven by: ", top_contributors$Dim2_Top, "\n",
  "Dim.3 is mainly driven by: ", top_contributors$Dim3_Top, "\n",
  "Dim.4 is mainly driven by: ", top_contributors$Dim4_Top, "\n",
  "Dim.5 is mainly driven by: ", top_contributors$Dim5_Top, "\n"
))

cat("\n--- BIOLOGICAL INTERPRETATION ---\n")
bio_interpretation <- function(params) {
  if (any(grepl("A|B", params))) {
    cat("- A and/or B dominate → photoperiod sensitivity (response to daylength).\n")
  }
  if (any(grepl("G|Th", params))) {
    cat("- G and/or Th dominate → thermal time / temperature sensitivity of development.\n")
  }
  if (any(grepl("Lc", params))) {
    cat("- Lc dominates → critical photoperiod threshold for reproductive transition.\n")
  }
  cat("\n")
}

cat("Dim.1 interpretation:\n")
bio_interpretation(top_contributors$Dim1_Top)
cat("Dim.2 interpretation:\n")
bio_interpretation(top_contributors$Dim2_Top)
cat("Dim.3 interpretation:\n")
bio_interpretation(top_contributors$Dim3_Top)

loadings_long <- melt(pca_loadings, id.vars = "Parameter",
                      variable.name = "Dimension",
                      value.name = "Loading")
loadings_plot <- loadings_long %>%
  filter(Dimension %in% c("Dim.1", "Dim.2"))

p_load <- ggplot(loadings_plot, aes(x = Parameter, y = Loading, fill = Dimension)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "PCA Loadings: Parameter Contributions (Dim.1 & Dim.2)",
       x = "Parameter", y = "Loading Value")
print(p_load)
save_plot(p_load, "PCA_Loadings_Dim1_Dim2")

### ============================================================
### STEP 7: Phenotypic Correlations (Overall + per Family)
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
    ggtitle("Phenotypic Correlations - All Families")
  print(ggp_all)
  save_plot(ggp_all, "Correlation_AllFamilies")
}

families <- unique(OptParams_complete$Family)
for(fam in families) {
  subset_fam <- OptParams_complete %>%
    filter(Family == fam) %>%
    select(all_of(numeric_params)) %>%
    drop_na()
  if (nrow(subset_fam) < 2) {
    cat("Skipping", fam, "- not enough taxa for correlation\n")
    next
  }
  ggp_fam <- ggpairs(subset_fam,
                     lower = list(continuous = wrap("points", alpha = 0.6)),
                     diag  = list(continuous = "densityDiag"),
                     upper = list(continuous = cor_sig_fn)) +
    ggtitle(paste("Phenotypic Correlations - Family:", fam))
  print(ggp_fam)
  save_plot(ggp_fam, paste0("Correlation_", fam))
}

### ============================================================
### STEP 7b: Family-wise G×A and B×Lc Correlations (Table)
### ============================================================

compute_family_correlation <- function(data, var1, var2) {
  cor_test <- tryCatch(cor.test(data[[var1]], data[[var2]], use = "complete.obs", method = "pearson"),
                       error = function(e) NULL)
  if (!is.null(cor_test)) {
    r_val <- cor_test$estimate
    p_val <- cor_test$p.value
  } else {
    r_val <- NA
    p_val <- NA
  }
  
  strength <- case_when(
    abs(r_val) >= 0.9 ~ "Very strong",
    abs(r_val) >= 0.7 ~ "Strong",
    abs(r_val) >= 0.5 ~ "Moderate",
    abs(r_val) >= 0.3 ~ "Weak",
    TRUE ~ "Negligible"
  )
  
  direction <- ifelse(r_val > 0, "Positive", "Negative")
  
  significance <- case_when(
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

for (fam in families) {
  fam_data <- OptParams_complete %>% filter(Family == fam)
  if (nrow(fam_data) > 3) {
    res_GA <- compute_family_correlation(fam_data, "G", "A") %>%
      mutate(Family = fam, n_Taxa = nrow(fam_data), Pair = "G × A")
    res_BLc <- compute_family_correlation(fam_data, "B", "Lc") %>%
      mutate(Family = fam, n_Taxa = nrow(fam_data), Pair = "B × Lc")
    results_GA <- bind_rows(results_GA, res_GA)
    results_BLc <- bind_rows(results_BLc, res_BLc)
  }
}

correlation_summary <- bind_rows(results_GA, results_BLc) %>%
  select(Family, Pair, n_Taxa, r, p, Strength, Direction, Significance, Interpretation) %>%
  arrange(Pair, desc(abs(r)))

cat("\n--- FAMILY-WISE G×A and B×Lc CORRELATIONS ---\n")
print(correlation_summary)

### ============================================================
### STEP 8: Genetic Correlation (Proxy via Family Means)
### ============================================================

param_cols <- c("G", "Th", "Lc", "A", "B")

family_means <- OptParams_complete %>%
  group_by(Family) %>%
  summarise(across(all_of(param_cols), mean, na.rm = TRUE))

genetic_cor <- cor(family_means[, param_cols], use = "pairwise.complete.obs", method = "pearson")
cor_test_results <- Hmisc::rcorr(as.matrix(family_means[, param_cols]), type = "pearson")
pvals <- cor_test_results$P

p_genetic <- ggcorrplot(
  genetic_cor,
  p.mat = pvals,
  hc.order = TRUE,
  type = "lower",
  lab = TRUE,
  lab_size = 4,
  title = "Genetic Correlation (Proxy) Based on Family Means",
  outline.color = "white",
  colors = c("blue", "white", "red")
)
print(p_genetic)
save_plot(p_genetic, "Genetic_Correlation")

### ============================================================
### STEP 8b: Within-Family (Residual) Correlations
### ============================================================

within_family_cor_list <- list()

for (fam in unique(OptParams_complete$Family)) {
  subset_data <- OptParams_complete %>%
    filter(Family == fam) %>%
    select(all_of(param_cols)) %>%
    drop_na()
  if (nrow(subset_data) > 3) {
    cor_res <- cor(subset_data, use = "pairwise.complete.obs", method = "pearson")
    within_family_cor_list[[fam]] <- cor_res
  }
}

within_family_summary <- do.call(rbind, lapply(names(within_family_cor_list), function(fam) {
  cor_mat <- within_family_cor_list[[fam]]
  data.frame(
    Family = fam,
    expand.grid(Var1 = rownames(cor_mat), Var2 = colnames(cor_mat)),
    Corr = as.vector(cor_mat)
  )
})) %>%
  filter(Var1 != Var2)

focus_pairs <- within_family_summary %>%
  filter((Var1 == "G" & Var2 == "A") | (Var1 == "B" & Var2 == "Lc"))

cat("\n--- WITHIN-FAMILY CORRELATIONS (G×A and B×Lc) ---\n")
print(focus_pairs)

avg_within_cor <- Reduce("+", within_family_cor_list) / length(within_family_cor_list)

p_resid <- ggcorrplot(
  avg_within_cor,
  hc.order = TRUE,
  type = "lower",
  lab = TRUE,
  lab_size = 4,
  title = "Average Within-Family Correlation (Residual Component)",
  outline.color = "white",
  colors = c("blue", "white", "red")
)
print(p_resid)
save_plot(p_resid, "Residual_Correlation")

### ============================================================
### STEP 8c: Compare Correlation Types (Phenotypic vs Genetic vs Residual)
### ============================================================

numeric_params <- c("G", "Th", "Lc", "A", "B")

phenotypic_corr <- cor(OptParams_complete[, numeric_params], use = "pairwise.complete.obs")
r_GA_pheno  <- phenotypic_corr["G", "A"]
r_BLc_pheno <- phenotypic_corr["B", "Lc"]

p_GA_pheno  <- cor.test(OptParams_complete$G,  OptParams_complete$A)$p.value
p_BLc_pheno <- cor.test(OptParams_complete$B,  OptParams_complete$Lc)$p.value

family_means_full <- OptParams_complete %>%
  group_by(Family) %>%
  summarise(across(all_of(numeric_params), mean, na.rm = TRUE)) %>%
  ungroup()

genetic_corr_full <- cor(family_means_full[, numeric_params], use = "pairwise.complete.obs")
r_GA_genetic  <- genetic_corr_full["G", "A"]
r_BLc_genetic <- genetic_corr_full["B", "Lc"]

p_GA_genetic  <- cor.test(family_means_full$G,  family_means_full$A)$p.value
p_BLc_genetic <- cor.test(family_means_full$B,  family_means_full$Lc)$p.value

within_corr_list <- list()
p_within_list    <- list()

for (fam in unique(OptParams_complete$Family)) {
  fam_data <- OptParams_complete %>%
    filter(Family == fam) %>%
    select(all_of(numeric_params)) %>%
    drop_na()
  
  if (nrow(fam_data) > 2) {
    cmat <- cor(fam_data, use = "pairwise.complete.obs")
    within_corr_list[[fam]] <- c(cmat["G", "A"], cmat["B", "Lc"])
    
    p_GA <- cor.test(fam_data$G, fam_data$A)$p.value
    p_BLc <- cor.test(fam_data$B, fam_data$Lc)$p.value
    p_within_list[[fam]] <- c(p_GA, p_BLc)
  }
}

within_corr_df <- do.call(rbind, within_corr_list)
p_within_df    <- do.call(rbind, p_within_list)
colnames(within_corr_df) <- c("GxA","BxLc")
colnames(p_within_df)    <- c("GxA","BxLc")

r_GA_resid  <- mean(within_corr_df[, "GxA"], na.rm = TRUE)
r_BLc_resid <- mean(within_corr_df[, "BxLc"], na.rm = TRUE)

p_GA_resid  <- mean(p_within_df[, "GxA"], na.rm = TRUE)
p_BLc_resid <- mean(p_within_df[, "BxLc"], na.rm = TRUE)

corr_summary <- data.frame(
  Pair        = c("GxA","BxLc"),
  Phenotypic  = c(r_GA_pheno,  r_BLc_pheno),
  Genetic     = c(r_GA_genetic, r_BLc_genetic),
  Residual    = c(r_GA_resid,  r_BLc_resid),
  p_Pheno     = c(p_GA_pheno,  p_BLc_pheno),
  p_Genetic   = c(p_GA_genetic, p_BLc_genetic),
  p_Residual  = c(p_GA_resid,  p_BLc_resid)
)

r_long <- corr_summary %>%
  select(Pair, Phenotypic, Genetic, Residual) %>%
  tidyr::pivot_longer(cols = c(Phenotypic, Genetic, Residual),
                      names_to = "Correlation_Type",
                      values_to = "r_value")

p_long <- corr_summary %>%
  select(Pair, p_Pheno, p_Genetic, p_Residual) %>%
  tidyr::pivot_longer(cols = c(p_Pheno, p_Genetic, p_Residual),
                      names_to = "p_name",
                      values_to = "p_value") %>%
  mutate(Correlation_Type = dplyr::recode(p_name,
                                          "p_Pheno"   = "Phenotypic",
                                          "p_Genetic" = "Genetic",
                                          "p_Residual"= "Residual")) %>%
  select(Pair, Correlation_Type, p_value)

corr_long <- dplyr::left_join(r_long, p_long,
                              by = c("Pair","Correlation_Type")) %>%
  mutate(sig = dplyr::case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01  ~ "**",
    p_value < 0.05  ~ "*",
    TRUE ~ ""
  ),
  label = paste0(sprintf("%.2f", r_value), sig))

corr_type_plot <- ggplot(corr_long, aes(x = Pair, y = r_value, fill = Correlation_Type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black") +
  geom_text(aes(label = label),
            position = position_dodge(width = 0.8),
            vjust = ifelse(corr_long$r_value > 0, -0.6, 1.3),
            size = 4.5) +
  scale_fill_manual(values = c("Phenotypic"="#6CA0DC",
                               "Genetic"="#E15759",
                               "Residual"="#59A14F")) +
  labs(title = "Comparison of Correlation Types (Phenotypic vs Genetic vs Residual)",
       x = "Parameter Pair", y = "Correlation Coefficient (r)",
       fill = "Correlation Type") +
  geom_hline(yintercept = 0, linetype="dashed", color="gray40") +
  theme_minimal(base_size = 14) +
  theme(legend.position="top",
        panel.grid.minor = element_blank(),
        plot.title = element_text(face="bold", size=16))
print(corr_type_plot)
save_plot(corr_type_plot, "Correlation_Type_Comparison")

### ============================================================
### STEP 9: Relate to Performance (Observed DTH)
### ============================================================

dth_vars <- c("DTH1", "DTH2", "DTH3")

family_perf <- OptParams_complete %>%
  group_by(Family) %>%
  summarise(across(c(G, Th, Lc, A, B, all_of(dth_vars)), mean, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(DTH_mean = rowMeans(across(all_of(dth_vars)), na.rm = TRUE))

param_model <- lm(DTH_mean ~ G + Th + Lc + A + B, data = family_perf)

cat("\n--- Regression Summary: Family Mean DTH ~ Parameters ---\n")
print(summary(param_model))

r2_param <- summary(param_model)$r.squared
cat("\nR² (Proportion of variance explained by parameters):", round(r2_param, 3), "\n")

param_importance <- broom::tidy(param_model) %>%
  filter(term != "(Intercept)") %>%
  mutate(Abs_Estimate = abs(estimate)) %>%
  arrange(desc(Abs_Estimate)) %>%
  rename(Parameter = term, Estimate = estimate, Std_Error = std.error,
         T_value = statistic, P_value = p.value)

cat("\n--- Parameter Importance (Sorted by Effect Size) ---\n")
print(param_importance)

pca_scores <- as.data.frame(pca_res$ind$coord)
pca_scores$Family <- pca_data$Family

family_pca <- pca_scores %>%
  group_by(Family) %>%
  summarise(across(starts_with("Dim"), mean, na.rm = TRUE))

pca_perf <- left_join(family_pca, family_perf, by = "Family")

pca_model <- lm(DTH_mean ~ Dim.1 + Dim.2 + Dim.3, data = pca_perf)

cat("\n--- Regression Summary: Family Mean DTH ~ PCA Components ---\n")
print(summary(pca_model))

r2_pca <- summary(pca_model)$r.squared
cat("\nR² (Proportion of variance explained by PCA components):", round(r2_pca, 3), "\n")

# --- Correlation of DTH_mean with parameters (family means) ---
cor_data <- family_perf[, c("DTH_mean", "G", "Th", "Lc", "A", "B")]
cor_res <- Hmisc::rcorr(as.matrix(cor_data), type = "pearson")
r_values <- cor_res$r
p_values <- cor_res$P

sig_stars <- ifelse(p_values < 0.001, "***",
                    ifelse(p_values < 0.01, "**",
                           ifelse(p_values < 0.05, "*", "")))

cor_table <- melt(r_values, varnames = c("Var1", "Var2"), value.name = "r")
p_table   <- melt(p_values, varnames = c("Var1", "Var2"), value.name = "p")
sig_table <- melt(sig_stars, varnames = c("Var1", "Var2"), value.name = "Significance")

cor_summary2 <- cor_table %>%
  left_join(p_table, by = c("Var1", "Var2")) %>%
  left_join(sig_table, by = c("Var1", "Var2")) %>%
  mutate(r = round(r, 3),
         p = signif(p, 3),
         Label = paste0(r, Significance))

cat("\n--- Correlation of DTH_mean with Parameters ---\n")
print(cor_summary2 %>% filter(Var1 == "DTH_mean" & Var2 != "DTH_mean"))

r_long2 <- cor_summary2 %>%
  filter(Var1 != Var2)

p_corr_heat <- ggplot(r_long2, aes(x = Var1, y = Var2, fill = r)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Label), size = 4, color = "black") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                       name = "r") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Correlation Matrix (r) with Significance",
       subtitle = "Stars indicate significance levels (* <0.05, ** <0.01, *** <0.001)",
       x = "", y = "")
print(p_corr_heat)
save_plot(p_corr_heat, "DTH_Parameter_Correlation_Heatmap")

param_cols <- c("G", "Th", "Lc", "A", "B")
for (pname in param_cols) {
  p_plot <- ggplot(family_perf, aes_string(x = pname, y = "DTH_mean", color = "Family")) +
    geom_point(size = 3) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
    theme_minimal(base_size = 14) +
    labs(title = paste("Relationship between", pname, "and Mean DTH"),
         x = pname, y = "Mean Observed DTH") +
    theme(legend.position = "none")
  print(p_plot)
  save_plot(p_plot, paste0("DTH_vs_", pname))
}

p_pca_dth <- ggplot(pca_perf, aes(x = Dim.1, y = DTH_mean, color = Family)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(title = "Relationship Between PC1 and Mean DTH",
       x = "PC1 (Multivariate Parameter Axis)",
       y = "Mean Observed DTH")
print(p_pca_dth)
save_plot(p_pca_dth, "PCA1_vs_DTH")

cat("\n--- Partial Regression Plots: saving as base-graphics PNG ---\n")
png(file.path(output_dir, "Partial_Regression_Plots.png"),
    width = 8, height = 6, units = "in", res = 300)
car::avPlots(param_model, ask = FALSE,
             main = "Partial Regression Plots: DTH_mean ~ G + Th + Lc + A + B")
dev.off()
cat("✅ Saved:", file.path(output_dir, "Partial_Regression_Plots.png"), "\n")

cat("\nINTERPRETATION GUIDE:\n")
cat("• Parameters with larger absolute coefficients (Estimate) have stronger effects on DTH.\n")
cat("• Positive estimate → increases DTH (later heading).\n")
cat("• Negative estimate → decreases DTH (earlier heading).\n")
cat("• High R² indicates parameters (or PCs) explain large proportion of heading variation.\n")
cat("• Partial regression plots show each parameter’s independent effect while controlling for others.\n")
cat("• PCA–DTH regression identifies whether combined multivariate variation in parameters predicts performance.\n")

### ============================================================
### EXTRA: 2A. Violin + Boxplots per Parameter by Family
### ============================================================

#for (pname in c("G", "Th", "Lc", "A", "B")) {
#  g_vio <- ggplot(OptParams_complete, aes(x = Family, y = .data[[pname]], fill = Family)) +
#    geom_violin(trim = FALSE, alpha = 0.6) +
#    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.5) +
#    theme_minimal(base_size = 13) +
#    theme(axis.text.x = element_text(angle = 45, hjust = 1),
#          legend.position = "none") +
#    labs(title = paste("Distribution of", pname, "by Family"),
#         y = pname, x = "Family")
#  print(g_vio)
#  save_plot(g_vio, paste0("ViolinBox_", pname))
#}

### ============================================================
### EXTRA: Flag Families with Unusual Within-Family Variance
### ============================================================

within_var_summary <- OptParams_complete %>%
  group_by(Family) %>%
  summarise(across(c(G, Th, Lc, A, B), var, na.rm = TRUE)) %>%
  pivot_longer(-Family, names_to = "Parameter", values_to = "Variance") %>%
  group_by(Parameter) %>%
  mutate(Zscore = scale(Variance))

outlier_families <- within_var_summary %>%
  filter(abs(Zscore) > 2)

cat("\n--- Families with unusually high/low within-family variance ---\n")
print(outlier_families)

g_var <- ggplot(within_var_summary, aes(x = Family, y = Variance, fill = Parameter)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(data = outlier_families,
            aes(label = "⚠", y = Variance + 0.05 * max(Variance)),
            position = position_dodge(width = 0.9),
            color = "red", size = 6) +
  theme_minimal(base_size = 13) +
  labs(title = "Within-Family Variance per Parameter (⚠ = outlier family)",
       y = "Variance", x = "Family") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(g_var)
save_plot(g_var, "Within_Family_Variance")

### ============================================================
### DONE
### ============================================================

cat("\n🎉 All figures have been saved to:\n", normalizePath(output_dir), "\n")
