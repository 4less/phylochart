library(phyloseq)
library(dplyr)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)

usage <- function() {
  cat(
    paste0(
      "Usage: Rscript export_ps_alpha_stats.R --input <phyloseq.rds> --output <alpha_stats.tsv> [--filter-treatment|--no-filter-treatment]\n",
      "Defaults: --input data/sylph_phyloseq.rds --output data/sylph_alpha_stats.tsv --filter-treatment\n"
    )
  )
  quit(status = 1)
}

opts <- list(
  input = "data/sylph_phyloseq.rds",
  output = "data/sylph_alpha_stats.tsv",
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

rank_candidates <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
ranks <- colnames(tax_table(ps)) %>% intersect(rank_candidates)
if (!length(ranks)) {
  stop("No supported taxonomic ranks found in tax_table.")
}

compute_alpha <- function(ps_obj) {
  otu <- phyloseq::otu_table(ps_obj)
  mat <- as(otu, "matrix")
  if (phyloseq::taxa_are_rows(otu)) {
    mat <- t(mat)
  }
  if (is.null(rownames(mat))) {
    stop("Sample names missing from OTU table.")
  }

  row_sums <- rowSums(mat)
  proportions <- sweep(mat, 1, row_sums, FUN = "/")

  observed <- rowSums(mat > 0)
  shannon <- apply(proportions, 1, function(p) {
    p <- p[is.finite(p) & p > 0]
    if (!length(p)) return(NA_real_)
    -sum(p * log(p))
  })
  simpson <- apply(proportions, 1, function(p) {
    p <- p[is.finite(p) & p > 0]
    if (!length(p)) return(NA_real_)
    1 - sum(p^2)
  })
  evenness <- ifelse(observed > 0, shannon / log(observed), NA_real_)

  data.frame(
    Sample = rownames(mat),
    Evenness = evenness,
    Shannon = shannon,
    Simpson = simpson,
    row.names = NULL
  )
}

alpha <- ranks %>%
  lapply(function(rank) {
    ps_rank <- tax_glom(ps, taxrank = rank, NArm = TRUE)
    metrics <- compute_alpha(ps_rank)
    metrics$Rank <- rank
    metrics
  }) %>%
  dplyr::bind_rows()

meta <- sample_data(ps) %>%
  data.frame() %>%
  tibble::rownames_to_column("Sample")

metrics <- alpha %>%
  select(Sample, Rank, Evenness, Shannon, Simpson) %>%
  pivot_longer(
    cols = c("Evenness", "Shannon", "Simpson"),
    names_to = "Type",
    values_to = "Value"
  )

alpha_long <- metrics %>%
  left_join(meta, by = "Sample")

meta_columns <- setdiff(colnames(meta), "Sample")
meta_columns <- meta_columns[meta_columns != "All"]

pairwise_results <- list()

pairwise_from_matrix <- function(pmat) {
  pairs <- which(!is.na(pmat), arr.ind = TRUE)
  if (!nrow(pairs)) return(data.frame())
  data.frame(
    Group1 = rownames(pmat)[pairs[, 1]],
    Group2 = colnames(pmat)[pairs[, 2]],
    PValue = pmat[pairs],
    stringsAsFactors = FALSE
  )
}

for (meta_col in meta_columns) {
  message("Testing metadata column: ", meta_col)
  column_values <- alpha_long[[meta_col]]
  if (all(is.na(column_values))) next
  alpha_long[[meta_col]] <- as.character(column_values)
  non_missing <- alpha_long %>%
    filter(!is.na(.data[[meta_col]])) %>%
    pull(.data[[meta_col]])
  if (length(non_missing) > 0 && length(unique(non_missing)) == length(non_missing)) {
    next
  }

  grouped <- alpha_long %>% group_by(Rank, Type)
  groups <- group_split(grouped)
  keys <- group_keys(grouped)

  results <- lapply(seq_along(groups), function(i) {
    df <- groups[[i]]
    key <- keys[i, ]
    df <- df %>%
      filter(!is.na(.data[[meta_col]]), !is.na(Value)) %>%
      mutate(Group = .data[[meta_col]])
    counts <- table(df$Group)
    keep_groups <- names(counts[counts >= 3])
    if (length(keep_groups) < 2) return(data.frame())
    df <- df %>% filter(.data$Group %in% keep_groups)
    test <- pairwise.wilcox.test(
      x = df$Value,
      g = df$Group,
      p.adjust.method = "none",
      exact = FALSE
    )
    pmat <- test$p.value
    res <- pairwise_from_matrix(pmat)
    if (!nrow(res)) return(data.frame())
    res$Rank <- key$Rank
    res$Type <- key$Type
    res$MetaColumn <- meta_col
    res$N1 <- as.integer(counts[res$Group1])
    res$N2 <- as.integer(counts[res$Group2])
    res$PAdj <- p.adjust(res$PValue, method = "BH")
    res$Test <- "pairwise_wilcoxon"
    res$PAdjustMethod <- "BH"
    res
  })

  pairwise_results[[meta_col]] <- bind_rows(results)
}

output <- bind_rows(pairwise_results)

readr::write_tsv(output, file = opts$output, quote = "none")
