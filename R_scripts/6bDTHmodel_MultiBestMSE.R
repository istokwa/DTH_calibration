# ============================================================
# 🔍 Extract taxa with multiple coarse exact minima (full details)
# ============================================================

library(readxl)
library(dplyr)
library(writexl)

# 1️⃣ Load the Coarse_Exact_Minima sheet
coarse_exact <- read_excel("results/GridSearch_G70_AllTaxa.xlsx",
                           sheet = "Coarse_Exact_Minima")

# 2️⃣ Count how many exact minima each taxon has
taxa_counts <- coarse_exact %>%
  group_by(Taxa) %>%
  summarise(ExactMinima_Count = n())

# 3️⃣ Identify taxa with multiple minima (≥2)
multi_minima_taxa <- taxa_counts %>%
  filter(ExactMinima_Count > 1)

# 4️⃣ Extract full rows for those taxa (with parameter values)
multi_minima_details <- coarse_exact %>%
  filter(Taxa %in% multi_minima_taxa$Taxa) %>%
  arrange(Taxa)

# 5️⃣ Print summary to console
cat("\n✅ Taxa with MULTIPLE exact minima in Coarse_Exact_Minima:\n")
print(multi_minima_taxa)
cat("\nTotal taxa with multiple exact minima:", nrow(multi_minima_taxa), "\n")

cat("\n✅ Corresponding parameter sets:\n")
print(multi_minima_details %>%
        select(Taxa, G, Th, Lc, A, B, MSE, starts_with("DTH"), starts_with("MODEL")))

# 6️⃣ Save results to Excel
write_xlsx(multi_minima_details, "results/Coarse_ExactMinima_MultiTaxaG70.xlsx")

cat("\n💾 Saved multiple-minima details to 'Coarse_ExactMinima_MultiTaxa.xlsx'\n")
