library(phyloseq)
library(dplyr)
library(vegan)

args <- commandArgs(trailingOnly = TRUE)

usage <- function() {
  cat(
    paste0(
      "Usage: Rscript export_ps_beta_stats.R --input <phyloseq.rds> --output <beta_stats.tsv> [--filter-treatment|--no-filter-treatment]\n",
      "Defaults: --input data/sylph_phyloseq.rds --output data/sylph_beta_stats.tsv --filter-treatment\n"
    )
  )
  quit(status = 1)
}

opts <- list(
  input = "data/sylph_phyloseq.rds",
  output = "data/sylph_beta_stats.tsv",
  filter_treatment = TRUE
)

idx <- 1
while (idx <= length(args)) {
  arg <- args[[idx]]
  if (arg == "--input") {
    idx <- idx + 1
    if (idx > length(args)) usage()
    opts$input <- args[[idx]]
  } else if (arg == "--output") {
    idx <- idx + 1
    if (idx > length(args)) usage()
    opts$output <- args[[idx]]
  } else if (arg == "--filter-treatment") {
    opts$filter_treatment <- TRUE
  } else if (arg == "--no-filter-treatment") {
    opts$filter_treatment <- FALSE
  } else if (arg %in% c("-h", "--help")) {
    usage()
  } else {
    cat("Unknown argument:", arg, "\n")
    usage()
  }
  idx <- idx + 1
}

ps <- readRDS(opts$input)
if (opts$filter_treatment) {
  meta_cols <- colnames(sample_data(ps))
  if ("Treatment" %in% meta_cols) {
    ps <- ps %>% subset_samples(!is.na(Treatment))
  } else {
    message("Warning: Treatment column not found; skipping treatment filter.")
  }
}

meta <- sample_data(ps) %>%
  data.frame() %>%
  tibble::rownames_to_column("Sample")

meta_columns <- setdiff(colnames(meta), "Sample")
meta_columns <- meta_columns[meta_columns != "All"]

methods <- c("bray", "jaccard")
results <- list()

for (method in methods) {
  message("Computing PERMANOVA: ", method)
  dist <- phyloseq::distance(ps, method = method)
  dist_mat <- as.matrix(dist)

  for (meta_col in meta_columns) {
    column_values <- meta[[meta_col]]
    if (all(is.na(column_values))) next
    column_values <- as.character(column_values)
    non_missing <- column_values[!is.na(column_values)]
    if (length(non_missing) > 0 && length(unique(non_missing)) == length(non_missing)) {
      next
    }
    df <- meta %>%
      select(Sample, Group = all_of(meta_col)) %>%
      filter(!is.na(Group))
    counts <- table(df$Group)
    keep_groups <- names(counts[counts >= 3])
    if (length(keep_groups) < 2) next
    df <- df %>% filter(Group %in% keep_groups)
    samples <- df$Sample
    dist_sub <- as.dist(dist_mat[samples, samples, drop = FALSE])
    adonis <- adonis2(dist_sub ~ Group, data = df, permutations = 999)
    if (!nrow(adonis)) next
    row <- adonis[1, ]
    results[[length(results) + 1]] <- data.frame(
      Metric = method,
      MetaColumn = meta_col,
      Df = row$Df,
      R2 = row$R2,
      F = row$F,
      PValue = row$`Pr(>F)`,
      stringsAsFactors = FALSE
    )
  }
}

output <- bind_rows(results)
if (nrow(output)) {
  output$PAdj <- p.adjust(output$PValue, method = "BH")
}
readr::write_tsv(output, file = opts$output, quote = "none")
