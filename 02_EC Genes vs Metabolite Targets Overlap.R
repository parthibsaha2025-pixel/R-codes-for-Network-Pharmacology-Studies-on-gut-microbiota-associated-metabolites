# ============================================================
# EC Genes vs Metabolite Targets - Overlap & Venn Diagram
# ============================================================
# Fully standalone: loads all required inputs itself via readxl,
# rather than assuming objects already exist in the R environment.
# ============================================================

if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")
if (!requireNamespace("VennDiagram", quietly = TRUE)) install.packages("VennDiagram")

library(readxl)
library(openxlsx)
library(VennDiagram)


# ------------------------------------------------------------
# SECTION 1: SETTINGS
# ------------------------------------------------------------
master_gene_file <- "PATH TO THE FILE"
ec_file          <- "PATH TO THE FILE"

venn_path    <- "PATH TO THE FILE"
overlap_path <- "PATH TO THE FILE"


# ------------------------------------------------------------
# SECTION 2: LOAD & CLEAN THE METABOLITE TARGET MASTER LIST
# ------------------------------------------------------------
if (!file.exists(master_gene_file)) {
  stop("Master gene list not found at: ", master_gene_file,
       " -- run metabolite_target_pipeline.R first to generate it.")
}

master_df <- read_excel(master_gene_file)
cat("Master gene list file rows (raw):", nrow(master_df), "\n")

master_raw <- trimws(master_df$gene)
master_na_like <- is.na(master_raw) | master_raw == "" |
  grepl("^n[\\s\\.\\/]*a$", master_raw, ignore.case = TRUE, perl = TRUE)

master_gene_list <- sort(unique(master_raw[!master_na_like]))

cat("Master list entries removed as NA/blank:", sum(master_na_like), "\n")
cat("Unique metabolite target genes loaded:", length(master_gene_list), "\n")


# ------------------------------------------------------------
# SECTION 3: LOAD & CLEAN THE EC GENE LIST
# ------------------------------------------------------------
if (!file.exists(ec_file)) {
  stop("EC gene file not found at: ", ec_file)
}

ec_df <- read_excel(ec_file)
cat("\nEC file rows (raw):", nrow(ec_df), "\n")

if (!"Symbol" %in% names(ec_df)) {
  stop("Column 'Symbol' not found in EC file. Columns present: ",
       paste(names(ec_df), collapse = ", "))
}

ec_raw <- trimws(ec_df$Symbol)
ec_na_like <- is.na(ec_raw) | ec_raw == "" |
  grepl("^n[\\s\\.\\/]*a$", ec_raw, ignore.case = TRUE, perl = TRUE)

ec_gene_list <- sort(unique(ec_raw[!ec_na_like]))

cat("EC entries removed as NA/blank:", sum(ec_na_like), "\n")
cat("Unique EC genes:", length(ec_gene_list), "\n")


# ------------------------------------------------------------
# SECTION 4: OVERLAP
# ------------------------------------------------------------
overlap_genes     <- sort(intersect(ec_gene_list, master_gene_list))
ec_only_genes     <- sort(setdiff(ec_gene_list, master_gene_list))
target_only_genes <- sort(setdiff(master_gene_list, ec_gene_list))

cat("\n---- OVERLAP SUMMARY ----\n")
cat("EC genes total:", length(ec_gene_list), "\n")
cat("Metabolite target genes total:", length(master_gene_list), "\n")
cat("Overlap (in both):", length(overlap_genes), "\n")
cat("EC-only:", length(ec_only_genes), "\n")
cat("Target-only:", length(target_only_genes), "\n")

cat("\nSanity check 1 - EC total reconstructs:",
    length(overlap_genes) + length(ec_only_genes) == length(ec_gene_list), "\n")
cat("Sanity check 2 - Target total reconstructs:",
    length(overlap_genes) + length(target_only_genes) == length(master_gene_list), "\n")


# ------------------------------------------------------------
# SECTION 5: VENN DIAGRAM
# ------------------------------------------------------------
venn.diagram(
  x = list("EC Genes" = ec_gene_list, "Metabolite Targets" = master_gene_list),
  filename = venn_path,
  imagetype = "png",
  height = 2000, width = 2000, resolution = 300,
  fill = c("#4E79A7", "#F28E2B"),
  alpha = 0.6,
  cex = 1.5,
  cat.cex = 1.3,
  cat.pos = c(-20, 20),
  main = "EC Genes vs Metabolite Targets",
  main.cex = 1.6,
  disable.logging = TRUE
)

cat("\nVenn diagram saved to:", venn_path, "\n")


# ------------------------------------------------------------
# SECTION 6: SAVE OVERLAP RESULTS
# ------------------------------------------------------------
out_wb <- createWorkbook()
addWorksheet(out_wb, "Overlap")
addWorksheet(out_wb, "EC_only")
addWorksheet(out_wb, "Target_only")

writeData(out_wb, "Overlap", data.frame(gene = overlap_genes))
writeData(out_wb, "EC_only", data.frame(gene = ec_only_genes))
writeData(out_wb, "Target_only", data.frame(gene = target_only_genes))

saveWorkbook(out_wb, overlap_path, overwrite = TRUE)
cat("Overlap results saved to:", overlap_path, "\n")


# ------------------------------------------------------------
# SECTION 7: FINAL SANITY CHECKS
# ------------------------------------------------------------
cat("\n---- FINAL SANITY CHECKS ----\n")
cat("1. Venn PNG exists:", file.exists(venn_path), "\n")
cat("2. Overlap Excel exists:", file.exists(overlap_path), "\n")

check_df <- read_excel(overlap_path, sheet = "Overlap")
cat("3. Rows in saved Overlap sheet:", nrow(check_df),
    "(should equal", length(overlap_genes), ") ->",
    ifelse(nrow(check_df) == length(overlap_genes), "OK", "MISMATCH"), "\n")

cat("4. Sample overlap genes:\n")
print(head(overlap_genes, 15))

cat("\n==== DONE ====\n")

