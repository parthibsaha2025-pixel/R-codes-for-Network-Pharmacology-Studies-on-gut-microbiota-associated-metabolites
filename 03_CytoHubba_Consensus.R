library(readr)
library(dplyr)

# ---- Config ----
NODE_TABLE_PATH <- "PATH TO THE FILE"
CYTOHUBBA_PATH  <- "PATH TO THE FILE"
TOP_N           <- 30

# ---- Read and map gene symbols ----
node_table <- read_csv(NODE_TABLE_PATH, show_col_types = FALSE)
hubba      <- read_csv(CYTOHUBBA_PATH, show_col_types = FALSE)

gene_map <- node_table %>%
  dplyr::select(name, `display name`) %>%
  dplyr::rename(node_name = name, gene = `display name`)

hubba_named <- hubba %>% left_join(gene_map, by = "node_name")
cat("Rows in cytoHubba export:", nrow(hubba), "\n")
cat("Rows mapped to a gene symbol:", sum(!is.na(hubba_named$gene)), "\n")

# ---- Flag top-N membership per method, then count consensus ----
hubba_methods <- c("MCC","DMNC","MNC","Degree","EPC","BottleNeck","EcCentricity",
                   "Closeness","Radiality","Betweenness","Stress","ClusteringCoefficient")

for (m in hubba_methods) {
  flag_col <- paste0(m, "_top", TOP_N)
  hubba_named[[flag_col]] <- rank(-hubba_named[[m]], ties.method = "min") <= TOP_N
}

flag_cols <- paste0(hubba_methods, "_top", TOP_N)
hubba_named$cytohubba_consensus <- rowSums(hubba_named[, flag_cols], na.rm = TRUE)

cytohubba_summary <- hubba_named %>%
  dplyr::select(gene, node_name, all_of(hubba_methods), cytohubba_consensus) %>%
  arrange(desc(cytohubba_consensus))

write_csv(cytohubba_summary, "PATH TO FILE")

cat("\nTop 20 by consensus (out of", length(hubba_methods), "methods):\n")
print(head(cytohubba_summary[, c("gene","cytohubba_consensus")], 20))

sum(is.na(hubba$ClusteringCoefficient))
