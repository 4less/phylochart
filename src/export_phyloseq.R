library(phyloseq)
library(purrr)
library(dplyr)

phyloseq_path <- "data/single_phyloseq.rds"
ps <- readRDS(phyloseq_path) %>%
  subset_samples(!is.na(Treatment))


ranks <- tax_table(ps) %>% 
  colnames() %>% 
  intersect(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))

tax <- tax_table(ps)

# Keep only ranks UP TO species
tax <- tax[, colnames(tax) %in% c("Kingdom","Phylum","Class","Order",
                                  "Family","Genus","Species")]

tax_table(ps) <- tax


taxon_rel_df <- ranks %>%
  map(~ {
    taxrank <- .x
    ps %>%
      transform_sample_counts(function(x) x / sum(x)) %>%
      tax_glom(taxrank = taxrank) %>%
      psmelt() %>%
      mutate(Taxon = .data[[taxrank]], Type = "RelativeAbundance", Rank = taxrank) %>%
      select(Sample, Taxon, Rank, Abundance, Type, meta_columns)
  }) %>%
  list_rbind()


taxon_abs_df <- ranks %>%
  map(~ {
    taxrank <- .x
    ps %>%
      tax_glom(taxrank = taxrank) %>%
      psmelt() %>%
      mutate(Taxon = .data[[taxrank]], Type = "AbsoluteAbundance", Rank = taxrank) %>%
      select(Sample, Taxon, Rank, Abundance, Type, meta_columns)
  }) %>%
  list_rbind()

readr::write_tsv(rbind(taxon_rel_df, taxon_abs_df), file = "data/taxa.tsv", quote="none")

ps %>%
  transform_sample_counts(function(x) x / sum(x)) %>%
  tax_glom(taxrank = "Genus") %>%
  psmelt() %>%
  mutate(Taxon = .data[["Genus"]], Type = "RelativeAbundance", Rank = "Genus") %>%
  select(Sample, Taxon, Rank, Abundance, Type, meta_columns) %>%
  filter(IDGroup == )
