### ============================================================
### STEP 0: Load Required Libraries
### ============================================================

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(car)       # for Levene's test
library(stats)     # for ANOVA, Kruskal-Wallis
library(DescTools) # optional descriptive stats
library(FSA)       # for Dunn’s post-hoc test
library(multcomp)  # for Tukey HSD post-hoc
library(rstatix)   # tidy post-hoc results formatting
library(agricolae)  # for CLD letters
library(ggpubr)     # for easier ggplot enhancements
library(lme4)      # for variance partitioning

setwd("D:/NU/Repository/GP_DTH-Rice/")

### ============================================================
### STEP 1: Data QC + Family-based Summary for Optimized Parameters
### ============================================================

#------------------------------------------------
# 1. Load data
#------------------------------------------------
OptParams <- read_excel("results/DVRparams.xlsx")

#------------------------------------------------
# 2. Define family mappings
#------------------------------------------------
family_map <- c(
  "WNAM_02" = "WNAM02",
  "WNAM_29" = "WNAM29",
  "WNAM_31" = "WNAM31",
  "WNAM_35" = "WNAM35",
  "WNAM_39" = "WNAM39",
  "WNAM_72" = "WNAM72",
  "WNAM_73" = "WNAM73"
)

#------------------------------------------------
# 3. Add Family column based on Taxa
#------------------------------------------------
OptParams <- OptParams %>%
  mutate(Family = substr(Taxa, 1, 7),
         Family = ifelse(Family %in% names(family_map),
                         family_map[Family],
                         Taxa))  # unique taxa get their own family label

#------------------------------------------------
# 4. Define numeric columns
#------------------------------------------------
numeric_cols <- c("G","Th","Lc","A","B","DVIStar","MSE",
                  "DTH1","DTH2","DTH3","MODEL1","MODEL2","MODEL3")

#------------------------------------------------
# 5. Data QC: Check for missing values
#------------------------------------------------
cat("\n--- MISSING VALUE CHECK ---\n")
na_summary <- colSums(is.na(OptParams[numeric_cols]))
print(na_summary)

if (any(na_summary > 0)) {
  cat("\nTaxa with missing data:\n")
  print(OptParams$Taxa[rowSums(is.na(OptParams[numeric_cols])) > 0])
} else {
  cat("✅ No missing values found.\n")
}

# Filter complete records only
OptParams_complete <- OptParams %>%
  filter(if_all(all_of(numeric_cols), ~ !is.na(.)))

#------------------------------------------------
# 6. Detect and flag outliers based on MSE
#------------------------------------------------
mse_threshold <- quantile(OptParams_complete$MSE, 0.95, na.rm = TRUE)
OptParams_complete <- OptParams_complete %>%
  mutate(FitQuality = ifelse(MSE > mse_threshold, "BadFit", "GoodFit"))

cat("\nMSE 95th percentile threshold:", round(mse_threshold, 3), "\n")
cat("Flagged bad fits:\n")
print(OptParams_complete$Taxa[OptParams_complete$FitQuality == "BadFit"])

# Count number of BadFit taxa
badfit_count <- OptParams_complete %>%
  filter(FitQuality == "BadFit") %>%
  summarise(Count = n())

cat("\nNumber of taxa flagged as BadFit:", badfit_count$Count, "\n")
cat("\nTaxa flagged as BadFit:\n")
print(OptParams_complete$Taxa[OptParams_complete$FitQuality == "BadFit"])

# Optional: Exclude bad fits
# OptParams_complete <- subset(OptParams_complete, FitQuality == "GoodFit")

#------------------------------------------------
# 7. Summary statistics per family (mean ± SD)
#------------------------------------------------
family_summary <- OptParams_complete %>%
  group_by(Family) %>%
  summarise(n = n(),
            across(all_of(numeric_cols),
                   ~ paste0(round(mean(.x, na.rm = TRUE), 2),
                            " ± ",
                            round(sd(.x, na.rm = TRUE), 2))),
            .groups = "drop")

#------------------------------------------------
# 8. Add overall summary
#------------------------------------------------
overall_summary <- OptParams_complete %>%
  summarise(n = n(),
            across(all_of(numeric_cols),
                   ~ paste0(round(mean(.x, na.rm = TRUE), 2),
                            " ± ",
                            round(sd(.x, na.rm = TRUE), 2)))) %>%
  mutate(Family = "Overall_NoGroup") %>%
  select(Family, everything())

#------------------------------------------------
# 9. Combine summaries
#------------------------------------------------
summary_table <- bind_rows(family_summary, overall_summary)
cat("\n--- FAMILY-LEVEL SUMMARY TABLE ---\n")
print(summary_table)

#------------------------------------------------
# 10. Quick visualization per parameter
#------------------------------------------------
#for (param in numeric_cols) {
#  if (is.numeric(OptParams_complete[[param]])) {
#    p <- ggplot(OptParams_complete, aes(x = Family, y = .data[[param]], fill = Family)) +
#      geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 1.5) +
#      theme_minimal(base_size = 13) +
#      theme(legend.position = "none",
#            axis.text.x = element_text(angle = 45, hjust = 1)) +
#      labs(title = paste("Parameter:", param),
#           x = "Family",
#           y = param)
#    print(p)
#  }

#}

cat("\n✅ QC and summary complete.\n")


### ============================================================
### STEP 2: Testing for Significant Differences Among Families
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
  
  if (length(unique(y)) <= 1) {
    param_test_results <- rbind(param_test_results,
                                data.frame(Parameter = param, Test_Type = NA, Statistic = NA,
                                           p_value = NA, Interpretation = "All values identical; test skipped"))
    next
  }
  
  # Normality and homogeneity
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
                              data.frame(Parameter = param, Test_Type = test_type,
                                         Statistic = stat, p_value = pval, Interpretation = interpretation))
}

cat("\n--- FAMILY COMPARISON TEST RESULTS ---\n")
print(param_test_results)


### ============================================================
### STEP 3: Testing for Significant Differences Within Families
### ============================================================

#------------------------------------------------
# Initialize storage
#------------------------------------------------
within_family_results <- data.frame(
  Family = character(),
  Parameter = character(),
  Test_Type = character(),
  Statistic = numeric(),
  p_value = numeric(),
  Interpretation = character(),
  stringsAsFactors = FALSE
)

#------------------------------------------------
# Loop through families and parameters
#------------------------------------------------
families <- unique(OptParams_complete$Family)

for (fam in families) {
  subset_data <- OptParams_complete %>% filter(Family == fam)
  
  for (param in numeric_cols) {
    y <- subset_data[[param]]
    taxa <- subset_data$Taxa
    
    # Skip if not enough variation
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
    
    # Test assumptions
    shapiro_p <- tryCatch({ shapiro.test(y)$p.value }, error = function(e) NA)
    levene_p <- tryCatch({ car::leveneTest(y ~ as.factor(taxa))$`Pr(>F)`[1] }, error = function(e) NA)
    
    parametric_ok <- !is.na(shapiro_p) && !is.na(levene_p) &&
      shapiro_p > 0.05 && levene_p > 0.05
    
    # Run appropriate test
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
### STEP 4: Post-Hoc Tests for Significant Parameters
### ============================================================

posthoc_results <- data.frame(
  Parameter = character(),
  Comparison = character(),
  p_value = numeric(),
  Adjusted_P = numeric(),
  Method = character(),
  stringsAsFactors = FALSE
)

for (i in 1:nrow(param_test_results)) {
  param <- param_test_results$Parameter[i]
  test_type <- param_test_results$Test_Type[i]
  pval <- param_test_results$p_value[i]
  
  if (is.na(pval) || pval >= 0.05) next
  
  cat(paste0("\nPerforming post-hoc for ", param, " (", test_type, ")...\n"))
  
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

# Optional: filter significant comparisons
significant_posthoc <- posthoc_results %>%
  filter(Adjusted_P < 0.05)
cat("\n--- SIGNIFICANT POST-HOC RESULTS ONLY ---\n")
print(significant_posthoc)

### ============================================================
### STEP 5: Diagnostic Plots with Post-Hoc Significance
### ============================================================

library(multcomp)
library(multcompView)

for (param in numeric_cols) {
  
  # Skip parameters with identical values
  if (length(unique(OptParams_complete[[param]])) <= 1) next
  
  y <- OptParams_complete[[param]]
  group <- as.factor(OptParams_complete$Family)
  test_type <- param_test_results$Test_Type[param_test_results$Parameter == param]
  
  cld_df <- NULL
  
  if (!is.na(test_type)) {
    
    if (test_type == "ANOVA") {
      # Tukey HSD via multcomp
      aov_model <- aov(y ~ group)
      tukey_res <- glht(aov_model, linfct = mcp(group = "Tukey"))
      cld_letters <- cld(tukey_res)
      cld_df <- data.frame(Family = names(cld_letters$mcletters$Letters),
                           CLD = cld_letters$mcletters$Letters,
                           stringsAsFactors = FALSE)
      
    } else if (test_type == "Kruskal-Wallis") {
      # Dunn test for pairwise comparisons
      dunn_res <- dunnTest(y ~ group, method = "bonferroni")$res
      
      # Build compact letter display using multcompView
      fams <- unique(group)
      p_matrix <- matrix(1, nrow = length(fams), ncol = length(fams))
      rownames(p_matrix) <- colnames(p_matrix) <- fams
      
      for (i in 1:nrow(dunn_res)) {
        comps <- strsplit(dunn_res$Comparison[i], " - ")[[1]]
        g1 <- comps[1]; g2 <- comps[2]
        p_matrix[g1, g2] <- dunn_res$P.adj[i]
        p_matrix[g2, g1] <- dunn_res$P.adj[i]
      }
      
      cld_letters <- multcompLetters(p_matrix, threshold = 0.05)
      cld_df <- data.frame(Family = names(cld_letters$Letters),
                           CLD = cld_letters$Letters,
                           stringsAsFactors = FALSE)
    }
  }
  
  # Prepare plot
  max_y <- max(y, na.rm = TRUE)
  y_range <- max_y - min(y, na.rm = TRUE)
  
  p <- ggplot(OptParams_complete, aes(x = Family, y = .data[[param]], fill = Family)) +
    geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 1.5) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste("Parameter:", param),
         x = "Family",
         y = param)
  
  # Add CLD letters above boxes
  if (!is.null(cld_df)) {
    cld_df <- cld_df %>% mutate(ypos = max_y + 0.1 * y_range)
    p <- p + geom_text(data = cld_df,
                       aes(x = Family, y = ypos, label = CLD),
                       inherit.aes = FALSE,
                       vjust = 0)
  }
  
  print(p)
}

cat("\n✅ Diagnostic plots with compact letter display complete.\n")


### ============================================================
### STEP 6: Phenotypic Correlations Among Parameters (Colored + Significance)
### ============================================================

library(GGally)
library(dplyr)
library(ggplot2)
library(Hmisc)  # for rcorr

#------------------------------------------------
# 1. Select numeric parameters
#------------------------------------------------
numeric_params <- c("G", "Th", "Lc", "A", "B")

# Custom function for upper triangle: correlation + significance
cor_sig_fn <- function(data, mapping, ...) {
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  
  # Compute correlation and p-value
  cor_res <- tryCatch(Hmisc::rcorr(x, y, type = "pearson"), error = function(e) NULL)
  
  if (!is.null(cor_res)) {
    r_val <- round(cor_res$r[1,2], 2)
    p_val <- cor_res$P[1,2]
    
    # Determine significance stars
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
  } else {
    ggplot() + theme_void()
  }
}

#------------------------------------------------
# 2. Overall correlation plot (all taxa)
#------------------------------------------------
subset_all <- OptParams_complete %>%
  select(all_of(numeric_params)) %>%
  drop_na()

if(nrow(subset_all) >= 2) {
  cat("\nGenerating overall phenotypic correlation plot...\n")
  
  ggp_all <- ggpairs(subset_all,
                     lower = list(continuous = wrap("points", alpha = 0.6)),
                     diag  = list(continuous = "densityDiag"),
                     upper = list(continuous = cor_sig_fn)) +
    ggtitle("Phenotypic Correlations - All Families")
  
  print(ggp_all)
} else {
  cat("Not enough data for overall correlation plot.\n")
}

#------------------------------------------------
# 3. Correlation plots per family
#------------------------------------------------
families <- unique(OptParams_complete$Family)

for(fam in families) {
  subset_fam <- OptParams_complete %>%
    filter(Family == fam) %>%
    select(all_of(numeric_params)) %>%
    drop_na()
  
  if(nrow(subset_fam) < 2) {
    cat("Skipping", fam, "- not enough taxa for correlation\n")
    next
  }
  
  cat("Generating correlation plot for family:", fam, "\n")
  
  ggp_fam <- ggpairs(subset_fam,
                     lower = list(continuous = wrap("points", alpha = 0.6)),
                     diag  = list(continuous = "densityDiag"),
                     upper = list(continuous = cor_sig_fn)) +
    ggtitle(paste("Phenotypic Correlations - Family:", fam))
  
  print(ggp_fam)
}


### ============================================================
### STEP 6b: Family-wise Parameter Correlations + Significance
### ============================================================

library(dplyr)

# Function to compute correlation with p-value and interpretation
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

# Prepare empty data frames for storing results
results_GA <- data.frame()
results_BLc <- data.frame()

# Get family list
families <- unique(OptParams_complete$Family)

# Loop over families
for (fam in families) {
  fam_data <- OptParams_complete %>% filter(Family == fam)
  
  if (nrow(fam_data) > 3) {
    # G vs A
    res_GA <- compute_family_correlation(fam_data, "G", "A") %>%
      mutate(Family = fam, n_Taxa = nrow(fam_data), Pair = "G × A")
    
    # B vs Lc
    res_BLc <- compute_family_correlation(fam_data, "B", "Lc") %>%
      mutate(Family = fam, n_Taxa = nrow(fam_data), Pair = "B × Lc")
    
    results_GA <- bind_rows(results_GA, res_GA)
    results_BLc <- bind_rows(results_BLc, res_BLc)
  }
}

# Combine both correlations into one table
correlation_summary <- bind_rows(results_GA, results_BLc) %>%
  select(Family, Pair, n_Taxa, r, p, Strength, Direction, Significance, Interpretation) %>%
  arrange(Pair, desc(abs(r)))

### ============================================================
### STEP 6: Multivariate Patterns — PCA and Clustering
### ============================================================

library(ggplot2)
library(dplyr)
library(FactoMineR)
library(factoextra)
library(cluster)

#------------------------------------------------
# 1. Select key numeric parameters for PCA
#------------------------------------------------
pca_data <- OptParams_complete %>%
  select(Taxa, Family, G, Th, Lc, A, B) %>%
  drop_na()

# Scale the numeric parameters (important for PCA)
pca_scaled <- scale(pca_data[, c("G", "Th", "Lc", "A", "B")])

#------------------------------------------------
# 2. Perform PCA
#------------------------------------------------
pca_res <- PCA(pca_scaled, graph = FALSE)

#------------------------------------------------
# 3. Visualize PCA results
#------------------------------------------------
fviz_pca_ind(
  pca_res,
  geom.ind = "point",
  col.ind = pca_data$Family,
  addEllipses = TRUE,
  palette = "Dark2",
  legend.title = "Family",
  repel = TRUE,
  title = "PCA of Optimized Parameters (G, Th, Lc, A, B)"
)

fviz_pca_var(
  pca_res,
  col.var = "contrib",
  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
  title = "Contribution of Parameters to Principal Components"
)

#------------------------------------------------
# 4. Hierarchical clustering by family
#------------------------------------------------
library(tidyverse)

family_means <- pca_data %>%
  group_by(Family) %>%
  summarise(across(c(G, Th, Lc, A, B), mean, na.rm = TRUE)) %>%
  column_to_rownames("Family")

dist_matrix <- dist(scale(family_means))
hc <- hclust(dist_matrix, method = "ward.D2")

plot(hc, main = "Hierarchical Clustering of Families (Ward.D2)",
     xlab = "Family", sub = "", cex = 1.2)

#------------------------------------------------
# 5. PCA Variance Summary
#------------------------------------------------
pca_summary <- data.frame(
  PC = paste0("PC", 1:length(pca_res$eig[,1])),
  Eigenvalue = round(pca_res$eig[,1], 3),
  Variance_percent = round(pca_res$eig[,2], 2),
  Cumulative_percent = round(pca_res$eig[,3], 2)
)

cat("\n--- PCA VARIANCE SUMMARY ---\n")
print(pca_summary)

#------------------------------------------------
# 6. Automatic PCA Interpretation
#------------------------------------------------
loadings <- as.data.frame(pca_res$var$coord)
colnames(loadings) <- paste0("PC", 1:ncol(loadings))
loadings$Variable <- rownames(loadings)

# Identify top contributing parameters for first 2 PCs
pc1_top <- loadings %>%
  arrange(desc(abs(PC1))) %>%
  slice(1:3) %>%
  pull(Variable)

pc2_top <- loadings %>%
  arrange(desc(abs(PC2))) %>%
  slice(1:3) %>%
  pull(Variable)

cat("\n--- PCA INTERPRETATION ---\n")
cat("PC1 (", round(pca_res$eig[1,2], 2), "% variance) is mainly driven by: ",
    paste(pc1_top, collapse = ", "), "\n", sep = "")
cat("PC2 (", round(pca_res$eig[2,2], 2), "% variance) is mainly driven by: ",
    paste(pc2_top, collapse = ", "), "\n", sep = "")

# Optional: create a short summary interpretation text
pca_text <- paste0(
  "PC1 explains ", round(pca_res$eig[1,2], 2), "% of total variance, mainly influenced by ",
  paste(pc1_top, collapse = ", "), 
  ". PC2 explains ", round(pca_res$eig[2,2], 2), "%, mainly influenced by ",
  paste(pc2_top, collapse = ", "), "."
)

cat("\nInterpretation summary:\n", pca_text, "\n")

### ============================================================
### STEP 6b: PCA Loadings + Biological Interpretation
### ============================================================

library(ggplot2)
library(reshape2)
library(dplyr)
library(factoextra)

# --- Extract loadings (parameter coordinates)
pca_loadings <- as.data.frame(pca_res$var$coord)
pca_loadings$Parameter <- rownames(pca_loadings)

cat("\n--- PCA LOADINGS TABLE ---\n")
# Round only numeric columns for viewing
print(round(pca_loadings[, sapply(pca_loadings, is.numeric)], 3))
cat("\nParameters:\n")
print(pca_loadings$Parameter)

# --- Rank parameters by contribution to each Dim
loading_strength <- pca_loadings %>%
  mutate(across(starts_with("Dim."), abs)) %>%
  arrange(desc(Dim.1)) %>%
  select(Parameter, Dim.1, Dim.2, Dim.3, Dim.4, Dim.5)

# Round only numeric columns for printing
cat("\n--- PARAMETERS CONTRIBUTING STRONGLY TO EACH DIMENSION ---\n")
loading_strength_rounded <- loading_strength
loading_strength_rounded[, sapply(loading_strength_rounded, is.numeric)] <-
  round(loading_strength_rounded[, sapply(loading_strength_rounded, is.numeric)], 3)
print(loading_strength_rounded)

# --- Automatic interpretation summary ---
cat("\n--- PCA INTERPRETATION SUMMARY ---\n")

# Determine top parameters per dimension
top_contributors <- loading_strength %>%
  summarise(
    Dim1_Top = paste(head(Parameter[order(-Dim.1)], 3), collapse = ", "),
    Dim2_Top = paste(head(Parameter[order(-Dim.2)], 3), collapse = ", "),
    Dim3_Top = paste(head(Parameter[order(-Dim.3)], 3), collapse = ", "),
    Dim4_Top = paste(head(Parameter[order(-Dim.4)], 3), collapse = ", "),
    Dim5_Top = paste(head(Parameter[order(-Dim.5)], 3), collapse = ", ")
  )

cat(paste0(
  "Dim.1 is mainly driven by: ", top_contributors$Dim1_Top, "\n",
  "Dim.2 is mainly driven by: ", top_contributors$Dim2_Top, "\n",
  "Dim.3 is mainly driven by: ", top_contributors$Dim3_Top, "\n",
  "Dim.4 is mainly driven by: ", top_contributors$Dim4_Top, "\n",
  "Dim.5 is mainly driven by: ", top_contributors$Dim5_Top, "\n"
))

# --- Biological interpretation of dominant PCs ---
cat("\n--- BIOLOGICAL INTERPRETATION ---\n")

bio_interpretation <- function(params) {
  if (any(grepl("A|B", params))) {
    cat("- Parameters A and B dominate → These reflect *photoperiod sensitivity* (response to daylength).\n")
  }
  if (any(grepl("G|Th", params))) {
    cat("- Parameters G and Th dominate → These relate to *thermal time and base temperature*, controlling growth rate and development speed.\n")
  }
  if (any(grepl("Lc", params))) {
    cat("- Parameter Lc contributes strongly → This indicates *critical photoperiod threshold*, the daylength boundary for reproductive transition.\n")
  }
  cat("\n")
}

cat("Dim.1 interpretation:\n")
bio_interpretation(top_contributors$Dim1_Top)
cat("Dim.2 interpretation:\n")
bio_interpretation(top_contributors$Dim2_Top)
cat("Dim.3 interpretation:\n")
bio_interpretation(top_contributors$Dim3_Top)

# --- Reshape for plotting
loadings_long <- melt(pca_loadings, id.vars = "Parameter",
                      variable.name = "Dimension",
                      value.name = "Loading")

# Focus only on Dim.1 and Dim.2
loadings_plot <- loadings_long %>%
  filter(Dimension %in% c("Dim.1", "Dim.2"))

# --- Visualization: Bar plot
ggplot(loadings_plot, aes(x = Parameter, y = Loading, fill = Dimension)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "PCA Loadings: Parameter Contributions (FactoMineR)",
       subtitle = "Positive and negative loadings show parameter direction along PC axes",
       x = "Parameter", y = "Loading Value") +
  scale_fill_brewer(palette = "Dark2")

# --- Optional: PCA Biplot (parameters as vectors)
#fviz_pca_biplot(pca_res, repel = TRUE,
#                title = "PCA Biplot: Parameter Loadings on Dim.1 & Dim.2")

### ============================================================
### STEP 6c (Safe version): PCA Clustering by Family (Active Only)
### ============================================================

library(factoextra)
library(dplyr)
library(ggplot2)

cat("\n--- PCA FAMILY CLUSTERING (Active Individuals Only) ---\n")

# ------------------------------------------------------------
# 1️⃣ Extract only active individuals from PCA
# ------------------------------------------------------------
active_ind_names <- rownames(pca_res$ind$coord)  # active taxa only

# Match only those taxa in your data
pca_plot_data <- OptParams_complete %>%
  filter(Taxa %in% active_ind_names) %>%
  arrange(match(Taxa, active_ind_names))

# Confirm alignment
if (length(active_ind_names) != nrow(pca_plot_data)) {
  stop("Mismatch: PCA active individuals and Family data still not aligned.")
}

# ------------------------------------------------------------
# 2️⃣ Visualize PCA individuals (active only)
# ------------------------------------------------------------
fviz_pca_ind(
  pca_res,
  geom.ind = "point",
  habillage = as.factor(pca_plot_data$Family),  # exactly matched length
  addEllipses = TRUE,
  ellipse.type = "norm",
  repel = TRUE,
  legend.title = "Family",
  palette = "Dark2",
  title = "PCA Clustering of Families (Active Individuals Only)"
)

# ------------------------------------------------------------
# 3️⃣ Compute centroid positions per family
# ------------------------------------------------------------
family_centroids <- pca_ind %>%
  group_by(Family) %>%
  summarise(Dim1 = mean(Dim.1), Dim2 = mean(Dim.2))

# ------------------------------------------------------------
# 4️⃣ Hierarchical clustering on family centroids
# ------------------------------------------------------------
dist_matrix <- dist(family_centroids[, c("Dim1", "Dim2")])
hc <- hclust(dist_matrix, method = "ward.D2")

# Dendrogram
plot(
  hc,
  main = "Hierarchical Clustering of Families (Ward’s Method)",
  xlab = "Family",
  ylab = "Euclidean Distance (PC Space)",
  sub = ""
)

# ------------------------------------------------------------
# 5️⃣ Optional: k-means clustering in PCA space
# ------------------------------------------------------------
set.seed(123)
k_value <- 3  # you can adjust number of clusters
kmeans_res <- kmeans(pca_ind[, c("Dim.1", "Dim.2")], centers = k_value, nstart = 25)
pca_ind$Cluster <- as.factor(kmeans_res$cluster)

# Visualization of k-means groups
ggplot(pca_ind, aes(x = Dim.1, y = Dim.2, color = Cluster, shape = Family)) +
  geom_point(size = 2, alpha = 0.8) +
  stat_ellipse(aes(color = Cluster), linetype = 2, level = 0.95) +
  theme_minimal(base_size = 13) +
  labs(title = paste("K-means Clustering (k =", k_value, ") in PCA Space"),
       x = "Dim.1", y = "Dim.2") +
  scale_color_brewer(palette = "Dark2")

# ------------------------------------------------------------
# 6️⃣ Interpret clustering
# ------------------------------------------------------------
cat("\n--- CLUSTERING INTERPRETATION ---\n")
cat("Families located near each other in PCA space have similar combinations of parameters.\n")
cat("Families far apart (or in distinct clusters) differ in their DTH-related parameter patterns.\n")
cat("Ellipses in the PCA plot indicate 95% confidence zones — small ellipses = more uniform family behavior.\n")
cat("You can use this to infer which families share similar temperature or photoperiod sensitivities.\n")

