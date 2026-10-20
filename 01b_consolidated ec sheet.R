if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
if (!requireNamespace("writexl", quietly = TRUE)) install.packages("writexl")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")

library(readxl)
library(writexl)
library(dplyr)

# ---- Read sources ----
GeneCards <- read_excel("PATH TO FILE")
CTD       <- read_excel("PATH TO FILE")
OMIM      <- read_excel("PATH TO FILE")

# ---- Standardize GeneCards column names ----
colnames(GeneCards) <- c("Details", "#", "Symbol", "Name", "Type", "Relevance Score", "Knowledge")

# GeneCards' Relevance Score / Knowledge are numeric; CTD's Reference Count
# column reads in as character (readxl fell back to text for at least one
# non-numeric cell), so GeneCards needs to match that type before bind_rows.
GeneCards$`Relevance Score` <- as.character(GeneCards$`Relevance Score`)
GeneCards$Knowledge         <- as.character(GeneCards$Knowledge)

# ---- Build CTD and OMIM into the same 7-column shape, all character ----
CTD_new <- data.frame(
  Details = rep("", nrow(CTD)),
  `#` = 1:nrow(CTD),
  Symbol = CTD$`Gene Symbol`,
  Name = rep("", nrow(CTD)),
  Type = rep("", nrow(CTD)),
  `Relevance Score` = as.character(CTD$`Inference Score`),
  Knowledge = as.character(CTD$`Reference Count`)
)
colnames(CTD_new) <- c("Details", "#", "Symbol", "Name", "Type", "Relevance Score", "Knowledge")

OMIM_new <- data.frame(
  Details = rep("", nrow(OMIM)),
  `#` = 1:nrow(OMIM),
  Symbol = OMIM$`Gene/Locus`,
  Name = OMIM$`Gene/Locus name`,
  Type = rep("", nrow(OMIM)),
  `Relevance Score` = rep("", nrow(OMIM)),
  Knowledge = rep("", nrow(OMIM))
)
colnames(OMIM_new) <- c("Details", "#", "Symbol", "Name", "Type", "Relevance Score", "Knowledge")

# ---- Combine ----
Combined_DB <- bind_rows(GeneCards, CTD_new, OMIM_new)
cat("Combined rows:", nrow(Combined_DB), "\n")
print(head(Combined_DB))
print(tail(Combined_DB))

Combined_DB <- Combined_DB %>% filter(!is.na(Symbol) & trimws(Symbol) != "")
Combined_DB$`#` <- 1:nrow(Combined_DB)

# Drop only rows with a genuinely blank Symbol, instead of blindly deleting
# row 1 — safer, since GeneCards is typically sorted by relevance and a
# hardcoded index-1 deletion risks silently dropping your top GeneCards hit
# if row 1 ever turns out to be real data rather than a blank artifact.
n_blank_symbol <- sum(trimws(Combined_DB$Symbol) == "")
cat("Rows with blank Symbol being dropped:", n_blank_symbol, "\n")
print(Combined_DB[trimws(Combined_DB$Symbol) == "", ])  # inspect before trusting the drop

Combined_DB <- Combined_DB %>% filter(trimws(Symbol) != "")
Combined_DB$`#` <- 1:nrow(Combined_DB)

cat("Final combined gene rows:", nrow(Combined_DB), "\n")
print(head(Combined_DB))

Combined_DB[is.na(Combined_DB)] <- ""

write_xlsx(Combined_DB, "PATH TO FILE")
cat("PATH TO FILE")
