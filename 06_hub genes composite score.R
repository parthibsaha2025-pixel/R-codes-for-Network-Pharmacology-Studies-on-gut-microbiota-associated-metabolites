library(readr)
library(dplyr)

# ---- Config ----
CYTOHUBBA_PATH   <- "PATH TO FILE"
PATHWAY_PATH     <- "PATH TO FILE"
CONVERGENCE_PATH <- "PATH TO FILE"
OUT_PATH         <- "PATH TO FILE"

CANDIDATE_GENES <- c("AKT1","GOT2", "PTGS2", "TP53", "EGFR", "MAPK3", 
                     "MAPK1", "ESR1", "PIK3CA", "SRC", "KRAS", "MAP2K1")

CYTOHUBBA_MAX <- 12   # total number of cytoHubba methods
PATHWAY_MAX   <- 6    # total number of curated EC pathways

# ---- Load the three scored tables ----
cytohubba_summary <- read_csv(CYTOHUBBA_PATH, show_col_types = FALSE)
pathway_summary   <- read_csv(PATHWAY_PATH, show_col_types = FALSE)
convergence       <- read_csv(CONVERGENCE_PATH, show_col_types = FALSE)

# Normalize metabolite convergence against the TRUE dataset-wide max
# (whatever gene has the highest n_metabolites across all 1,238 genes,
# not just among the candidates below)
metabolite_max <- max(convergence$n_metabolites, na.rm = TRUE)
cat("Metabolite convergence normalized against dataset max:", metabolite_max, "\n")

# ---- Build the composite table ----
composite_table <- tibble(gene = CANDIDATE_GENES) %>%
  left_join(cytohubba_summary %>% dplyr::select(gene, cytohubba_consensus), by = "gene") %>%
  left_join(pathway_summary   %>% dplyr::select(gene, ec_pathway_count, ec_pathways), by = "gene") %>%
  left_join(convergence       %>% dplyr::select(gene, n_metabolites), by = "gene") %>%
  mutate(
    cytohubba_norm   = cytohubba_consensus / CYTOHUBBA_MAX,
    ec_pathway_norm  = ec_pathway_count / PATHWAY_MAX,
    metabolite_norm  = n_metabolites / metabolite_max,
    composite_score  = cytohubba_norm + ec_pathway_norm + metabolite_norm
  ) %>%
  arrange(desc(composite_score))

print(composite_table %>%
        dplyr::select(gene, cytohubba_consensus, ec_pathway_count, n_metabolites, composite_score))

write_csv(composite_table, OUT_PATH)
cat("\nSaved to:", OUT_PATH, "\n")




