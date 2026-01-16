### ============================================================
### CORRELATION ANALYSIS ONLY (GGPAIRS OUTPUT)
### ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(GGally)
  library(Hmisc)     # rcorr()
})

setwd("D:/NU/Repository/GP_DTH-Rice/")

experimentID <- "enviTest"
runDate <- "(12-24-2025)"

# --- Output directory ---
output_dir <- "results/3 - Correlations/"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

filename <- paste0(runDate, "-Correlation-", experimentID)

save_plot <- function(plot_obj, filename, width = 7, height = 6, dpi = 300) {
  filepath <- file.path(output_dir, paste0(filename, ".png"))
  ggsave(filepath, plot = plot_obj, width = width, height = height, dpi = dpi)
  cat("✅ Saved:", filepath, "\n")
}

# --- Load data (your optimized parameter sheet) ---
OptParams <- read_excel("data/processed/(12-23-2025)-enviTest.xlsx", sheet = "Finer_Best_Params")

# --- Family column (keep if you want per-family plots) ---
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
  mutate(
    Family = substr(Taxa, 1, 7),
    Family = ifelse(Family %in% names(family_map), family_map[Family], Family)
  )

# --- Correlation variables (match your figure) ---
numeric_params <- c("G", "Th", "Lc", "A", "B")

# --- Clean complete cases for correlation plotting ---
OptParams_complete <- OptParams %>%
  select(Family, all_of(numeric_params)) %>%
  drop_na()

# --- Upper panel function: r + significance stars (and color like your figure) ---
cor_sig_fn <- function(data, mapping, ...) {
  x <- GGally::eval_data_col(data, mapping$x)
  y <- GGally::eval_data_col(data, mapping$y)
  
  cor_res <- tryCatch(Hmisc::rcorr(x, y, type = "pearson"),
                      error = function(e) NULL)
  if (is.null(cor_res)) return(ggplot() + theme_void())
  
  r_val <- round(cor_res$r[1, 2], 2)
  p_val <- cor_res$P[1, 2]
  sig_star <- ifelse(p_val < 0.001, "***",
                     ifelse(p_val < 0.01,  "**",
                            ifelse(p_val < 0.05, "*", "")))
  
  txt_col <- ifelse(r_val > 0.5, "red",
                    ifelse(r_val < -0.5, "blue", "black"))
  
  ggplot(data = data, mapping = mapping) +
    annotate("text", x = 0.5, y = 0.5,
             label = paste0(r_val, sig_star),
             size = 6, color = txt_col) +
    theme_void()
}

### ============================================================
### 1) PHENOTYPIC CORRELATIONS — ALL FAMILIES (ONE FIGURE)
### ============================================================

subset_all <- OptParams_complete %>%
  select(all_of(numeric_params))

ggp_all <- ggpairs(
  subset_all,
  lower = list(continuous = wrap("points", alpha = 0.6, size = 1)),
  diag  = list(continuous = "densityDiag"),
  upper = list(continuous = cor_sig_fn)
) +
  ggtitle(paste0("Correlations of DVR Parameters - ", experimentID)) +
  theme(
    # ✅ parameter names (G, Th, Lc, A, B) live in facet strips:
    strip.text = element_text(size = 14, face = "bold"),
    # optional readability tweaks:
    axis.text  = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  )

print(ggp_all)
save_plot(ggp_all, filename, width = 7, height = 6)

cat("\n🎉 Saved correlation figures to:\n", normalizePath(output_dir), "\n")
