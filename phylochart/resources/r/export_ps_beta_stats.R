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

ranks <- tax_table(ps) %>%
  colnames() %>%
  intersect(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
if (!length(ranks)) {
  ranks <- "All"
}

methods <- c("bray", "jaccard", "clr+euclidean")
results <- list()

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
    message("Computing PERMANOVA: ", method, " (", rank, ")")
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
        Test = "global",
        Comparison = "All",
        Metric = method,
        Rank = rank,
        MetaColumn = meta_col,
        Df = row$Df,
        R2 = row$R2,
        F = row$F,
        PValue = row$`Pr(>F)`,
        stringsAsFactors = FALSE
      )

      groups <- sort(unique(as.character(df$Group)))
      if (length(groups) < 2) next
      combos <- combn(groups, 2, simplify = FALSE)
      for (pair in combos) {
        df_pair <- df %>% filter(Group %in% pair)
        if (length(unique(df_pair$Group)) < 2) next
        samples_pair <- df_pair$Sample
        dist_pair <- as.dist(dist_mat[samples_pair, samples_pair, drop = FALSE])
        adonis_pair <- adonis2(dist_pair ~ Group, data = df_pair, permutations = 999)
        if (!nrow(adonis_pair)) next
        pair_row <- adonis_pair[1, ]
        results[[length(results) + 1]] <- data.frame(
          Test = "pairwise",
          Comparison = paste(pair, collapse = " vs "),
          Metric = method,
          Rank = rank,
          MetaColumn = meta_col,
          Df = pair_row$Df,
          R2 = pair_row$R2,
          F = pair_row$F,
          PValue = pair_row$`Pr(>F)`,
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

output <- bind_rows(results)
if (nrow(output)) {
  output$PAdj <- p.adjust(output$PValue, method = "BH")
}
readr::write_tsv(output, file = opts$output, quote = "none")
