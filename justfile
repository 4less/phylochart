build:
  python scripts/build_phylochart.py data/taxa.tsv

build-sylph:
  python scripts/build_phylochart.py --taxa data/sylph_taxa.tsv --alpha data/sylph_alpha.tsv --alpha-stats data/sylph_alpha_stats.tsv --beta data/sylph_beta.tsv --beta-stats data/sylph_beta_stats.tsv

sylph-taxa:
  Rscript src/export_ps_taxa.R --input data/sylph_phyloseq.rds --output data/sylph_taxa.tsv

sylph-alpha:
  Rscript src/export_ps_alpha.R --input data/sylph_phyloseq.rds --output data/sylph_alpha.tsv

sylph-alpha-stats:
  Rscript src/export_ps_alpha_stats.R --input data/sylph_phyloseq.rds --output data/sylph_alpha_stats.tsv

sylph-beta:
  Rscript src/export_ps_beta.R --input data/sylph_phyloseq.rds --output data/sylph_beta.tsv

sylph-beta-stats:
  Rscript src/export_ps_beta_stats.R --input data/sylph_phyloseq.rds --output data/sylph_beta_stats.tsv
