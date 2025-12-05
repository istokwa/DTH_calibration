# ============================================================
# 🔍 Extract taxa with only one coarse exact minimum (full details)
# ============================================================

library(readxl)
library(dplyr)
library(writexl)

# 1️⃣ Load the Coarse_Exact_Minima sheet
coarse_exact <- read_excel("results/DVRparams-DVRparams-G_DTH2.xlsx",
                           sheet = "Finer_Exact_Minima")

# 2️⃣ Count and filter taxa with exactly one exact minimum
single_minima_taxa <- coarse_exact %>%
  group_by(Taxa) %>%
  summarise(ExactMinima_Count = n()) %>%
  filter(ExactMinima_Count == 1)

# 3️⃣ Extract the detailed rows for those taxa
single_minima_details <- coarse_exact %>%
  filter(Taxa %in% single_minima_taxa$Taxa) %>%
  arrange(Taxa)

# 4️⃣ Show summary in console
cat("\n✅ Taxa with only ONE exact minimum in Coarse_Exact_Minima:\n")
print(unique(single_minima_details$Taxa))
cat("\nTotal taxa with single exact minima:", 
    length(unique(single_minima_details$Taxa)), "\n")

# 5️⃣ Save results (only details)
write_xlsx(single_minima_details, "results/MultiTaxa-G_DTH2.xlsx")

cat("\n💾 Saved single-minima details to 'Coarse_ExactMinima_SingleTaxa.xlsx'\n")
