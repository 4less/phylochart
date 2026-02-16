build:
  phylochart build --taxa data/taxa.tsv --output phylochart.html

build-sylph:
  phylochart build --taxa data/sylph_taxa.tsv --alpha data/sylph_alpha.tsv --alpha-stats data/sylph_alpha_stats.tsv --beta data/sylph_beta.tsv --beta-stats data/sylph_beta_stats.tsv --pcoa data/sylph_pcoa.tsv --output phylochart.html

build-single:
  phylochart build --taxa data/single_taxa.tsv --alpha data/single_alpha.tsv --alpha-stats data/single_alpha_stats.tsv --beta data/single_beta.tsv --beta-stats data/single_beta_stats.tsv --pcoa data/single_pcoa.tsv --output phylochart.html

sylph-taxa:
  phylochart taxa --input data/sylph_phyloseq.rds --output data/sylph_taxa.tsv

sylph-alpha:
  phylochart alpha --input data/sylph_phyloseq.rds --output data/sylph_alpha.tsv --output-stats data/sylph_alpha_stats.tsv

sylph-alpha-stats:
  phylochart alpha-stats --input data/sylph_phyloseq.rds --output data/sylph_alpha_stats.tsv

sylph-beta:
  phylochart beta --input data/sylph_phyloseq.rds --output data/sylph_beta.tsv --output-stats data/sylph_beta_stats.tsv

sylph-beta-stats:
  phylochart beta-stats --input data/sylph_phyloseq.rds --output data/sylph_beta_stats.tsv

sylph-pcoa:
  phylochart pcoa --input data/sylph_phyloseq.rds --output data/sylph_pcoa.tsv


single-taxa:
  phylochart taxa --input data/single_phyloseq.rds --output data/single_taxa.tsv

single-alpha:
  phylochart alpha --input data/single_phyloseq.rds --output data/single_alpha.tsv --output-stats data/single_alpha_stats.tsv

single-alpha-stats:
  phylochart alpha-stats --input data/single_phyloseq.rds --output data/single_alpha_stats.tsv

single-beta:
  phylochart beta --input data/single_phyloseq.rds --output data/single_beta.tsv --output-stats data/single_beta_stats.tsv

single-beta-stats:
  phylochart beta-stats --input data/single_phyloseq.rds --output data/single_beta_stats.tsv

single-pcoa:
  phylochart pcoa --input data/single_phyloseq.rds --output data/single_pcoa.tsv


sylph-all:
  phylochart all --input data/sylph_phyloseq.rds --output-dir data/sylph_report

single-all:
  phylochart all --input data/single_phyloseq.rds --output-dir data/single_report
