build:
  python scripts/build_phylochart.py data/taxa.tsv

build-sylph:
  python scripts/build_phylochart.py --taxa data/sylph_taxa.tsv --alpha data/sylph_alpha.tsv --alpha-stats data/sylph_alpha_stats.tsv --beta data/sylph_beta.tsv --beta-stats data/sylph_beta_stats.tsv --pcoa data/sylph_pcoa.tsv

build-single:
  python scripts/build_phylochart.py --taxa data/single_taxa.tsv --alpha data/single_alpha.tsv --alpha-stats data/single_alpha_stats.tsv --beta data/single_beta.tsv --beta-stats data/single_beta_stats.tsv --pcoa data/single_pcoa.tsv

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

sylph-pcoa:
  Rscript src/export_ps_pcoa.R --input data/sylph_phyloseq.rds --output data/sylph_pcoa.tsv


single-taxa:
  Rscript src/export_ps_taxa.R --input data/single_phyloseq.rds --output data/single_taxa.tsv

single-alpha:
  Rscript src/export_ps_alpha.R --input data/single_phyloseq.rds --output data/single_alpha.tsv

single-alpha-stats:
  Rscript src/export_ps_alpha_stats.R --input data/single_phyloseq.rds --output data/single_alpha_stats.tsv

single-beta:
  Rscript src/export_ps_beta.R --input data/single_phyloseq.rds --output data/single_beta.tsv

single-beta-stats:
  Rscript src/export_ps_beta_stats.R --input data/single_phyloseq.rds --output data/single_beta_stats.tsv

single-pcoa:
  Rscript src/export_ps_pcoa.R --input data/single_phyloseq.rds --output data/single_pcoa.tsv


sylph-all:
  just sylph-taxa
  just sylph-alpha
  just sylph-alpha-stats
  just sylph-beta
  just sylph-beta-stats
  just sylph-pcoa

single-all:
  just single-taxa
  just single-alpha
  just single-alpha-stats
  just single-beta
  just single-beta-stats
  just single-pcoa
