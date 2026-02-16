library(phyloseq)
library(dplyr)
library(vegan)

args <- commandArgs(trailingOnly = TRUE)

usage <- function() {
  cat(
    paste0(
      "Usage: Rscript export_ps_pcoa.R --input <phyloseq.rds> --output <pcoa.tsv> [--filter-treatment|--no-filter-treatment]\n",
      "Defaults: --input data/sylph_phyloseq.rds --output data/sylph_pcoa.tsv --filter-treatment\n"
    )
  )
  quit(status = 1)
}

opts <- list(
  input = "data/sylph_phyloseq.rds",
  output = "data/sylph_pcoa.tsv",
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
    message("Computing ordination: ", method, " (", rank, ")")
    ps_rank <- build_rank_ps(rank)
    if (nsamples(ps_rank) < 3 || ntaxa(ps_rank) < 2) {
      message("Skipping ", method, " (", rank, "): insufficient samples or taxa.")
      next
    }
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
    if (!nrow(dist_mat)) next
    if (any(!is.finite(dist_mat)) || all(dist_mat == 0)) {
      message("Skipping ", method, " (", rank, "): distance matrix degenerate.")
      next
    }

    k <- min(3, nrow(dist_mat) - 1)
    if (k >= 2) {
      pcoa <- tryCatch(
        cmdscale(dist, k = k, eig = TRUE),
        error = function(e) NULL
      )
      if (is.null(pcoa) || is.null(pcoa$points)) {
        message("Skipping PCoA for ", method, " (", rank, "): cmdscale failed.")
        next
      }
      points <- as.data.frame(pcoa$points)
      if (!ncol(points) || ncol(points) < 2) {
        message("Skipping PCoA for ", method, " (", rank, "): no axes returned.")
        next
      }
      colnames(points) <- paste0("Axis", seq_len(ncol(points)))
      samples <- rownames(points)

      eig_vals <- pcoa$eig
      eig_pos <- eig_vals[eig_vals > 0]
      total <- sum(eig_pos)
      if (!is.finite(total) || total <= 0) {
        message("Skipping PCoA for ", method, " (", rank, "): no positive eigenvalues.")
        next
      }
      var_vals <- eig_vals / total

      coords <- data.frame(
        Sample = samples,
        Metric = method,
        Rank = rank,
        Method = "PCoA",
        Axis1 = points$Axis1,
        Axis2 = if ("Axis2" %in% names(points)) points$Axis2 else NA_real_,
        Axis3 = if ("Axis3" %in% names(points)) points$Axis3 else NA_real_,
        Var1 = var_vals[1],
        Var2 = var_vals[2],
        Var3 = var_vals[3],
        row.names = NULL,
        stringsAsFactors = FALSE
      )
      all_rows[[paste(method, rank, "PCoA", sep = "_")]] <- coords
    }

    if (nrow(dist_mat) >= 4) {
      nmds <- tryCatch(
        metaMDS(dist, k = 2, trymax = 50, autotransform = FALSE, trace = 0),
        error = function(e) NULL
      )
      if (!is.null(nmds)) {
        points <- as.data.frame(nmds$points)
        colnames(points) <- c("Axis1", "Axis2")
        samples <- rownames(points)
        coords <- data.frame(
          Sample = samples,
          Metric = method,
          Rank = rank,
          Method = "NMDS",
          Axis1 = points$Axis1,
          Axis2 = points$Axis2,
          Axis3 = NA_real_,
          Var1 = NA_real_,
          Var2 = NA_real_,
          Var3 = NA_real_,
          row.names = NULL,
          stringsAsFactors = FALSE
        )
        all_rows[[paste(method, rank, "NMDS", sep = "_")]] <- coords
      }
    }
  }
}

output <- bind_rows(all_rows)
readr::write_tsv(output, file = opts$output, quote = "none")
