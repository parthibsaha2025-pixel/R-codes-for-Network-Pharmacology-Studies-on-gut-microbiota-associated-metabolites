# One-time install if needed:
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "KEGGREST"))
library(readxl)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(KEGGREST)

# ---- Config ----
CORE_TARGETS_PATH <- "PATH TO THE FILE"
EC_PATHWAYS <- c(
  "Endometrial cancer", "PI3K-Akt signaling pathway", "Estrogen signaling pathway",
  "p53 signaling pathway", "Apoptosis", "MAPK signaling pathway"
)

# ---- Map gene symbols -> Entrez IDs ----
core_targets <- read_excel(CORE_TARGETS_PATH)
gene_symbols <- unique(trimws(core_targets$gene))
cat("Core target set size:", length(gene_symbols), "\n")

id_map <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
cat("Mapped to Entrez ID:", nrow(id_map), "out of", length(gene_symbols), "\n")
entrez_ids <- unique(id_map$ENTREZID)

# ---- GO (Biological Process) ----
go_bp <- enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, ont = "BP",
                  pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
cat("Significant GO:BP terms:", nrow(as.data.frame(go_bp)), "\n")

# ---- KEGG (still used for the dotplot / general enrichment overview) ----
kegg <- enrichKEGG(gene = entrez_ids, organism = "hsa",
                   pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
kegg_readable <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
kegg_df <- as.data.frame(kegg_readable)
cat("Significant KEGG pathways:", nrow(kegg_df), "\n")

write.csv(as.data.frame(go_bp), "PATH TO FILE", row.names = FALSE)
write.csv(kegg_df, "PATH TO FILE", row.names = FALSE)

png("PATH TO FILE", width = 2000, height = 1600, res = 200)
dotplot(go_bp, showCategory = 20) + ggplot2::ggtitle("GO Biological Process")
dev.off()

png("PATH TO FILE", width = 2000, height = 1600, res = 200)
dotplot(kegg_readable, showCategory = 20) + ggplot2::ggtitle("KEGG Pathways")
dev.off()

# ------------------------------------------------------------
# EC-pathway plausibility, using KEGG's FIXED pathway gene-set
# definitions (independent of this run's enrichment significance)
# ------------------------------------------------------------

# 1. Resolve each curated pathway name to its hsa pathway ID by looking
#    it up directly against KEGG's full human pathway list (avoids
#    hand-typed hsa IDs and any exact-string-match fragility)
all_hsa_pathways <- keggList("pathway", "hsa")
cleaned_pw_names <- sub(" - Homo sapiens.*$", "", all_hsa_pathways)

match_pathway_id <- function(target_name) {
  hit <- names(all_hsa_pathways)[cleaned_pw_names == target_name]
  if (length(hit) == 0) {
    hit <- names(all_hsa_pathways)[grepl(target_name, cleaned_pw_names, ignore.case = TRUE)]
  }
  hit
}

ec_pathway_ids <- sapply(EC_PATHWAYS, match_pathway_id)
cat("\nResolved pathway IDs (check every entry has exactly one hsa ID):\n")
print(ec_pathway_ids)

# 2. Pull the FULL gene membership for each pathway directly from KEGG
pathway_membership <- lapply(ec_pathway_ids, function(pid) {
  if (length(pid) == 0 || is.na(pid[1])) return(character(0))
  links <- keggLink("hsa", pid[1])
  entrez <- sub("^hsa:", "", links)
  syms <- suppressWarnings(mapIds(org.Hs.eg.db, keys = entrez,
                                  keytype = "ENTREZID", column = "SYMBOL"))
  unique(na.omit(syms))
})
names(pathway_membership) <- EC_PATHWAYS

cat("\nGene-set size per curated pathway (sanity check - none should be 0):\n")
print(sapply(pathway_membership, length))

# 3. Score every core target gene against these fixed gene sets
ec_pathway_count <- sapply(gene_symbols, function(g) {
  sum(sapply(pathway_membership, function(genes) g %in% genes))
})
ec_pathway_which <- sapply(gene_symbols, function(g) {
  hits <- names(pathway_membership)[sapply(pathway_membership, function(genes) g %in% genes)]
  paste(hits, collapse = "; ")
})

pathway_summary <- data.frame(
  gene = gene_symbols,
  ec_pathway_count = ec_pathway_count,
  ec_pathways = ec_pathway_which,
  stringsAsFactors = FALSE
) %>% arrange(desc(ec_pathway_count))

write.csv(pathway_summary, "PATH TO FILE", row.names = FALSE)

cat("\nTop 20 by EC-pathway plausibility:\n")
print(head(pathway_summary, 20))
