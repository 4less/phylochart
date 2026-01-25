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

ranks <- tax_table(ps) %>%
  colnames() %>%
  intersect(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
if (!length(ranks)) {
  ranks <- "All"
}

methods <- c("bray", "jaccard", "clr+euclidean")
all_rows <- list()

build_rank_ps <- function(rank) {
  if (rank == "All") {
    return(ps)
  }
  tax_glom(
    ps,
    taxrank = rank,
    NArm = TRUE,
    bad_empty = c(NA, "", " ", "\t")
  )
}

for (method in methods) {
  for (rank in ranks) {
    message("Computing beta distances: ", method, " (", rank, ")")
    ps_rank <- build_rank_ps(rank)
    if (method == "clr+euclidean") {
      ps_rank <- transform_sample_counts(ps_rank, function(x) {
        x <- x + 1
        logx <- log(x)
        logx - mean(logx)
      })
      dist <- phyloseq::distance(ps_rank, method = "euclidean")
    } else {
      dist <- phyloseq::distance(ps_rank, method = method)
    }
    mat <- as.matrix(dist)
    if (!nrow(mat)) next
    samples <- rownames(mat)
    idx_pairs <- which(upper.tri(mat), arr.ind = TRUE)
    rows <- data.frame(
      Sample1 = samples[idx_pairs[, 1]],
      Sample2 = samples[idx_pairs[, 2]],
      Metric = method,
      Rank = rank,
      Distance = mat[idx_pairs],
      row.names = NULL,
      stringsAsFactors = FALSE
    )
    all_rows[[paste(method, rank, sep = "_")]] <- rows
  }
}

output <- bind_rows(all_rows)
readr::write_tsv(output, file = opts$output, quote = "none")
