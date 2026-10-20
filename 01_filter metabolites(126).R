getwd()
# -----------------------------------------------------------------------
# Filter a folder of metabolite subfolders down to a target list
# -----------------------------------------------------------------------
# Edit these three paths/settings, then run the script.

metabolite_list_file <- "PATH TO THE FILE"   # your list of 126 metabolites
# (a .csv/.txt with one metabolite
# name per line, or a column
# named "metabolite")
source_dir <- "PATH TO THE FILE"
output_dir <- "PATH TO THE FILE"     # will be created; matching
# folders are COPIED here

case_insensitive <- TRUE   # match "Glucose" == "glucose"
trim_whitespace  <- TRUE   # strip leading/trailing spaces from names
strip_parens     <- TRUE   # "Lithocholic acid (LCA)" -> "Lithocholic acid"
fuzzy_fallback   <- TRUE   # suggest close matches for anything still unmatched
fuzzy_max_dist   <- 0.1    # max edit distance (fraction of chars) for suggestions

# -----------------------------------------------------------------------
# 1. Read the metabolite list
# -----------------------------------------------------------------------
ext <- tolower(tools::file_ext(metabolite_list_file))

if (ext %in% c("csv", "txt", "tsv")) {
  # Try reading as a table; if it has a header with a metabolite-like
  # column, use that column, otherwise treat the first column as the list.
  df <- read.csv(metabolite_list_file, header = TRUE, stringsAsFactors = FALSE,
                 sep = if (ext == "tsv") "\t" else ",")
  name_col <- grep("metabolite|name|compound", names(df), ignore.case = TRUE, value = TRUE)
  if (length(name_col) > 0) {
    metabolites <- df[[name_col[1]]]
  } else {
    metabolites <- df[[1]]
  }
} else if (ext %in% c("xlsx", "xls")) {
  if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
  df <- readxl::read_excel(metabolite_list_file)
  name_col <- grep("metabolite|name|compound", names(df), ignore.case = TRUE, value = TRUE)
  metabolites <- if (length(name_col) > 0) df[[name_col[1]]] else df[[1]]
} else {
  stop("Unsupported file type for metabolite list: ", ext)
}

metabolites <- as.character(metabolites)
metabolites <- metabolites[!is.na(metabolites) & metabolites != ""]
if (trim_whitespace) metabolites <- trimws(metabolites)

cat("Loaded", length(metabolites), "metabolites from list.\n")

# -----------------------------------------------------------------------
# 2. List the existing subfolders
# -----------------------------------------------------------------------
all_folders <- list.dirs(source_dir, full.names = FALSE, recursive = FALSE)
cat("Found", length(all_folders), "subfolders in source directory.\n")

clean_name <- function(x) {
  if (strip_parens) x <- gsub("\\s*\\([^)]*\\)", "", x)   # drop "(LCA)" etc.
  if (trim_whitespace) x <- trimws(x)
  if (case_insensitive) x <- tolower(x)
  x
}

folder_key <- clean_name(all_folders)
met_key    <- clean_name(metabolites)

# -----------------------------------------------------------------------
# 3. Match (exact match on cleaned names)
# -----------------------------------------------------------------------
matched_lgl  <- met_key %in% folder_key
match_idx    <- all_folders[folder_key %in% met_key[matched_lgl]]

cat("Matched", sum(matched_lgl), "of", length(metabolites), "metabolites to folders (exact, after cleaning).\n")

missing <- metabolites[!matched_lgl]

# -----------------------------------------------------------------------
# 3b. Fuzzy suggestions for anything still unmatched (review manually!)
# -----------------------------------------------------------------------
if (length(missing) > 0 && fuzzy_fallback) {
  cat("\nNo exact match for", length(missing), "metabolites. Nearest folder name(s):\n")
  for (i in seq_along(missing)) {
    mkey <- clean_name(missing[i])
    dists <- adist(mkey, folder_key)[1, ]
    max_d <- ceiling(fuzzy_max_dist * nchar(mkey))
    best  <- which(dists <= max(max_d, 1))
    if (length(best) > 0) {
      best <- best[order(dists[best])]
      cat(sprintf("  '%s'  ->  %s\n", missing[i],
                  paste0(all_folders[best], " (dist=", dists[best], ")", collapse = "; ")))
    } else {
      cat(sprintf("  '%s'  ->  no close match found\n", missing[i]))
    }
  }
  cat("\nReview the suggestions above. If any are correct, either rename the\n",
      "folder to match your list, or add the metabolite name to the list\n",
      "using the folder's exact spelling, then re-run this script.\n", sep = "")
} else if (length(missing) > 0) {
  cat("\nMetabolites with NO matching folder (", length(missing), "):\n", sep = "")
  print(missing)
}

extra <- setdiff(all_folders, match_idx)
cat("\n", length(extra), "folders will be excluded (not in your list).\n")

# -----------------------------------------------------------------------
# 4. Copy matching folders to output_dir
# -----------------------------------------------------------------------
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

for (folder in match_idx) {
  src  <- file.path(source_dir, folder)
  dest <- file.path(output_dir, folder)
  file.copy(src, output_dir, recursive = TRUE, overwrite = TRUE)
}

cat("\nDone. Copied", length(match_idx), "folders into:\n", normalizePath(output_dir), "\n")


# ------------------------------------------------------------
# SECTION 8: MASTER GENE LIST WITH SOURCE + PROBABILITY
#            (applying probability cutoffs: STP >= 0.5, STITCH >= 0.4)
# ------------------------------------------------------------
output_folder <- "PATH TO FILE"
# 8.1 - gather metabolite, gene, source, probability from every consolidated file
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
cat("Total raw entries:", nrow(raw_data), "\n")

# 8.2 - idempotent NA/"N/A"/blank gene cleaning, always recomputed from raw
clean_raw_data <- function(df) {
  g <- trimws(df$gene)
  na_like <- is.na(g) | g == "" | grepl("^n[\\s\\.\\/]*a$", g, ignore.case = TRUE, perl = TRUE)
  df2 <- df[!na_like, ]
  df2$gene <- trimws(df2$gene)
  df2
}

data_clean <- clean_raw_data(raw_data)
cat("Entries removed as NA/N.A./blank:", nrow(raw_data) - nrow(data_clean), "\n")

table(data_clean$source)
# 8.3 - apply probability cutoffs, per source
STP_CUTOFF    <- 0.5
STITCH_CUTOFF <- 0.4

passes_cutoff <- (data_clean$source == "STP"    & data_clean$probability >= STP_CUTOFF) |
  (data_clean$source == "STITCH" & data_clean$probability >= STITCH_CUTOFF)

n_before_cutoff <- nrow(data_clean)
n_genes_before  <- length(unique(data_clean$gene))

data_filtered <- data_clean[passes_cutoff, ]

cat("\n---- CUTOFF FILTERING ----\n")
cat("Rows before cutoff:", n_before_cutoff, "| after:", nrow(data_filtered),
    "| removed:", n_before_cutoff - nrow(data_filtered), "\n")
cat("  removed STP rows below", STP_CUTOFF, ":",
    sum(data_clean$source == "STP" & data_clean$probability < STP_CUTOFF, na.rm = TRUE), "\n")
cat("  removed STITCH rows below", STITCH_CUTOFF, ":",
    sum(data_clean$source == "STITCH" & data_clean$probability < STITCH_CUTOFF, na.rm = TRUE), "\n")

genes_dropped_entirely <- setdiff(unique(data_clean$gene), unique(data_filtered$gene))
cat("Genes with ALL occurrences below cutoff (dropped from list entirely):",
    length(genes_dropped_entirely), "\n")

# 8.4 - source summary per gene (across ALL metabolites, post-cutoff): STP / STITCH / Both
source_summary <- aggregate(
  source ~ gene, data = data_filtered,
  FUN = function(s) { u <- unique(s); if (length(u) > 1) "Both" else u }
)

# 8.5 - probability: within each (metabolite, gene) pair, if both STP and
# STITCH rows survive the cutoff, keep only the STP probability for that metabolite
dedup_by_metabolite_gene <- do.call(rbind, lapply(
  split(data_filtered, list(data_filtered$metabolite, data_filtered$gene), drop = TRUE),
  function(grp) {
    if ("STP" %in% grp$source) grp <- grp[grp$source == "STP", , drop = FALSE]
    grp[1, , drop = FALSE]
  }
))
cat("Rows after per-metabolite STP-priority dedup:", nrow(dedup_by_metabolite_gene), "\n")

# 8.6 - collapse across metabolites to one probability per gene (max = strongest prediction)
prob_summary <- aggregate(
  probability ~ gene, data = dedup_by_metabolite_gene,
  FUN = function(p) max(p, na.rm = TRUE)
)
names(prob_summary) <- c("gene", "max_probability")

# 8.7 - merge source + probability into the final master table
master_df <- merge(source_summary, prob_summary, by = "gene", all = TRUE)
master_df <- master_df[order(master_df$gene), ]

cat("\nUnique genes in final master list (post-cutoff):", nrow(master_df), "\n")
print(table(master_df$source))

# 8.8 - save
out_path <- file.path("PATH TO FILE")
write.xlsx(master_df, out_path)
cat("Saved to:", out_path, "\n")

# ------------------------------------------------------------
# SECTION 9: SANITY CHECKS
# ------------------------------------------------------------
cat("\n---- SANITY CHECKS ----\n")

# 9.1 - gene count vs pre-cutoff count (expect a decrease, not a fixed number anymore)
cat("1. Genes before cutoff:", n_genes_before,
    "| after cutoff:", nrow(master_df),
    "| dropped:", n_genes_before - nrow(master_df),
    "(should equal", length(genes_dropped_entirely), ") ->",
    ifelse(n_genes_before - nrow(master_df) == length(genes_dropped_entirely), "OK", "CHECK"), "\n")

# 9.2 - every surviving probability actually meets its source's cutoff
#       (checked against the per-metabolite-gene rows that fed into max_probability,
#        since max_probability itself is a max across sources/metabolites)
below_cutoff_survived <- sum(
  (dedup_by_metabolite_gene$source == "STP" & dedup_by_metabolite_gene$probability < STP_CUTOFF) |
    (dedup_by_metabolite_gene$source == "STITCH" & dedup_by_metabolite_gene$probability < STITCH_CUTOFF)
)
cat("2. Rows feeding into the master list that are below their cutoff:", below_cutoff_survived,
    "->", ifelse(below_cutoff_survived == 0, "OK", "CHECK"), "\n")

# 9.3 - every gene has a probability value, values sensibly bounded
n_missing_prob <- sum(is.na(master_df$max_probability))
out_of_range <- sum(master_df$max_probability < 0 | master_df$max_probability > 1, na.rm = TRUE)
cat("3. Missing probability values:", n_missing_prob,
    "| out of [0,1] range:", out_of_range,
    "->", ifelse(n_missing_prob == 0 && out_of_range == 0, "OK", "CHECK"), "\n")

# 9.4 - spot-check a few genes and their probabilities
sample_genes <- intersect(c("CA2", "CA1", "PPARG"), master_df$gene)
cat("4. Spot-check sample genes:\n")
print(master_df[master_df$gene %in% sample_genes, ])

# 9.5 - show a sample of genes dropped entirely (useful to sanity-check the cutoff wasn't too aggressive)
cat("5. Sample of genes dropped entirely (all occurrences below cutoff):\n")
print(head(sort(genes_dropped_entirely), 15))

# 9.6 - saved file matches in-memory result
file_ok <- file.exists(out_path)
written <- if (file_ok) read.xlsx(out_path) else NULL
cat("6. File exists:", file_ok,
    "| rows match:", !is.null(written) && nrow(written) == nrow(master_df), "\n")

summary(data_clean$probability[data_clean$source == "STITCH"])

