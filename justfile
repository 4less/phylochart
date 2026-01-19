build:
  python scripts/build_phylochart.py data/taxa.tsv

build-sylph:
  python scripts/build_phylochart.py data/sylph_taxa.tsv

sylph_taxa:
  Rscript src/export_phyloseq.R --input data/sylph_phyloseq.rds --output data/sylph_taxa.tsv
