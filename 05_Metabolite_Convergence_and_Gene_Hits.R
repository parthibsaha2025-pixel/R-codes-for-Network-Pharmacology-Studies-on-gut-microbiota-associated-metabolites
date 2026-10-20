library(openxlsx)
library(dplyr)

# ============================================================
# CONFIG
# ============================================================
output_folder <- "PATH TO FILE"
STP_CUTOFF    <- 0.5
STITCH_CUTOFF <- 0.4

# Genes you want a closer look at — used for both the convergence
# preview and the metabolite-hit lookup below
CANDIDATE_GENES <- c("AKT1","GOT2", "PTGS2", "TP53", "EGFR", "MAPK3", 
                     "MAPK1", "ESR1", "PIK3CA", "SRC", "KRAS", "MAP2K1")

# ============================================================
# 1. Re-gather metabolite, gene, source, probability from every consolidated file
# ============================================================
consolidated_files <- list.files(output_folder, pattern = "_Consolidated\\.xlsx$",
                                 recursive = TRUE, full.names = TRUE)
cat("Consolidated files found:", length(consolidated_files), "\n")

raw_data <- do.call(rbind, lapply(consolidated_files, function(f) {
  df <- read.xlsx(f)
  metabolite <- sub("_Consolidated\\.xlsx$", "", basename(f))
  data.frame(
    metabolite = metabolite,
    gene = df$gene,
    source = df$source,
    probability = as.numeric(df$probability),
    stringsAsFactors = FALSE
  )
}))

# ============================================================
# 2. Clean NA/blank genes
# ============================================================
g <- trimws(raw_data$gene)
na_like <- is.na(g) | g == "" | grepl("^n[\\s\\.\\/]*a$", g, ignore.case = TRUE, perl = TRUE)
data_clean <- raw_data[!na_like, ]
data_clean$gene <- trimws(data_clean$gene)

# ============================================================
# 3. Apply cutoffs
# ============================================================
passes_cutoff <- (data_clean$source == "STP"    & data_clean$probability >= STP_CUTOFF) |
  (data_clean$source == "STITCH" & data_clean$probability >= STITCH_CUTOFF)
data_filtered <- data_clean[passes_cutoff, ]
cat("Rows after cutoff:", nrow(data_filtered), "\n")

# ============================================================
# 4. Convergence: how many DISTINCT metabolites hit each gene
#    (a metabolite counts once even if it hit the gene via both STP and STITCH)
# ============================================================
convergence <- data_filtered %>%
  distinct(metabolite, gene) %>%
  count(gene, name = "n_metabolites") %>%
  arrange(desc(n_metabolites))
cat("Genes with convergence data:", nrow(convergence), "\n")

write.csv(convergence, "PATH TO FILE", row.names = FALSE)

cat("\nConvergence distribution:\n")
print(table(convergence$n_metabolites))

cat("\nTop 70 genes by convergence:\n")
print(head(convergence, 70))

cat("\nConvergence for candidate genes:\n")
print(convergence %>% filter(gene %in% CANDIDATE_GENES) %>% arrange(desc(n_metabolites)))

# ============================================================
# 5. Which specific metabolites hit the candidate genes
# ============================================================
gene_hits <- data_filtered %>%
  filter(gene %in% CANDIDATE_GENES) %>%
  arrange(gene, desc(probability))

cat("\nMetabolite-level hits for candidate genes:\n")
print(gene_hits)

write.csv(gene_hits, "PATH TO FILE", row.names = FALSE)

