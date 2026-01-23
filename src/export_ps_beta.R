library(phyloseq)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

usage <- function() {
  cat(
    paste0(
      "Usage: Rscript export_ps_beta.R --input <phyloseq.rds> --output <beta.tsv> [--filter-treatment|--no-filter-treatment]\n",
      "Defaults: --input data/sylph_phyloseq.rds --output data/sylph_beta.tsv --filter-treatment\n"
    )
  )
  quit(status = 1)
}

opts <- list(
  input = "data/sylph_phyloseq.rds",
  output = "data/sylph_beta.tsv",
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

methods <- c("bray", "jaccard")
all_rows <- list()

for (method in methods) {
  message("Computing beta distances: ", method)
  dist <- phyloseq::distance(ps, method = method)
  mat <- as.matrix(dist)
  if (!nrow(mat)) next
  samples <- rownames(mat)
  idx_pairs <- which(upper.tri(mat), arr.ind = TRUE)
  rows <- data.frame(
    Sample1 = samples[idx_pairs[, 1]],
    Sample2 = samples[idx_pairs[, 2]],
    Metric = method,
    Distance = mat[idx_pairs],
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  all_rows[[method]] <- rows
}

output <- bind_rows(all_rows)
readr::write_tsv(output, file = opts$output, quote = "none")
