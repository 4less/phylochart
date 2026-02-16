library(phyloseq)
library(purrr)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

usage <- function() {
  cat(
    paste0(
      "Usage: Rscript export_ps_taxa.R --input <phyloseq.rds> --output <taxa.tsv> [--filter-treatment|--no-filter-treatment] [--types relative|absolute|both] [--drop-zeros|--keep-zeros]\n",
      "Defaults: --input data/sylph_phyloseq.rds --output data/sylph_taxa.tsv --filter-treatment --types both --drop-zeros\n"
    )
  )
  quit(status = 1)
}

opts <- list(
  input = "data/sylph_phyloseq.rds",
  output = "data/sylph_taxa.tsv",
  filter_treatment = TRUE,
  types = "both",
  drop_zeros = TRUE
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
  } else if (arg == "--types") {
    idx <- idx + 1
    if (idx > length(args)) usage()
    opts$types <- args[[idx]]
  } else if (arg == "--drop-zeros") {
    opts$drop_zeros <- TRUE
  } else if (arg == "--keep-zeros") {
    opts$drop_zeros <- FALSE
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

tax <- tax_table(ps)

# Keep only ranks UP TO species
tax <- tax[, colnames(tax) %in% c("Kingdom","Phylum","Class","Order",
                                  "Family","Genus","Species")]

tax_table(ps) <- tax

meta_columns <- colnames(sample_data(ps))

allowed_types <- c("relative", "absolute", "both")
if (!(opts$types %in% allowed_types)) {
  cat("Invalid --types value:", opts$types, "\n")
  usage()
}

taxon_rel_df <- tibble()
taxon_abs_df <- tibble()


has_absolute <- {
  ma <- phyloseq::sample_sums(ps) %>% mean()
  ma > 110 & ma < 90
}

message(phyloseq::sample_sums(ps) %>% mean())

if (opts$types %in% c("relative", "both")) {
  taxon_rel_df <- ranks %>%
    map(~ {
      taxrank <- .x
      message("RelativeAbundance: ", taxrank)
      ps %>%
        transform_sample_counts(function(x) x / sum(x)) %>%
        tax_glom(taxrank = taxrank) %>%
        psmelt() %>%
        { if (opts$drop_zeros) filter(., Abundance > 0) else . } %>%
        mutate(Taxon = .data[[taxrank]], Type = "RelativeAbundance", Rank = taxrank) %>%
        select(Sample, Taxon, Rank, Abundance, Type, all_of(meta_columns))
    }) %>%
    list_rbind()
}

if (opts$types %in% c("absolute", "both") & has_absolute) {
  taxon_abs_df <- ranks %>%
    map(~ {
      taxrank <- .x
      message("AbsoluteAbundance: ", taxrank)
      ps %>%
        tax_glom(taxrank = taxrank) %>%
        psmelt() %>%
        { if (opts$drop_zeros) filter(., Abundance > 0) else . } %>%
        mutate(Taxon = .data[[taxrank]], Type = "AbsoluteAbundance", Rank = taxrank) %>%
        select(Sample, Taxon, Rank, Abundance, Type, all_of(meta_columns))
    }) %>%
    list_rbind()
}

readr::write_tsv(rbind(taxon_rel_df, taxon_abs_df), file = opts$output, quote = "none")




taxrank <- "Species"
test <- ps %>%
  transform_sample_counts(function(x) x / sum(x)) %>%
  tax_glom(taxrank=taxrank, NArm=TRUE, bad_empty=c(NA, "", " ", "\t")) %>%
  psmelt() %>%
  mutate(Taxon = .data[[taxrank]], Type = "RelativeAbundance", Rank = taxrank) %>%
  select(Sample, Taxon, Rank, Abundance, Type, all_of(meta_columns))

test %>% filter()

test %>% group_by(Sample) %>%
  summarize(total = sum(Abundance))


