# phylochart

Install from GitHub (single line):

```bash
python -m pip install "git+https://github.com/4less/phylochart.git"
```

Local editable install (optional):

```bash
python -m pip install -e .
```

CLI:

```bash
phylochart alpha --input phyloseq.rds --output alpha.tsv --output-stats alpha_stats.tsv
phylochart beta --input phyloseq.rds --output beta.tsv --output-stats beta_stats.tsv
phylochart pcoa --input phyloseq.rds --output pcoa.tsv
phylochart taxa --input phyloseq.rds --output taxa.tsv
phylochart build --taxa taxa.tsv --alpha alpha.tsv --alpha-stats alpha_stats.tsv --beta beta.tsv --beta-stats beta_stats.tsv --pcoa pcoa.tsv --output phylochart.html
phylochart all --input phyloseq.rds --output-dir out
```

`phylochart all` writes:
- `out/stats/taxa.tsv`
- `out/stats/alpha.tsv`
- `out/stats/alpha_stats.tsv`
- `out/stats/beta.tsv`
- `out/stats/beta_stats.tsv`
- `out/stats/pcoa.tsv`
- `out/phylochart.html`

Example (end-to-end):

```bash
phylochart all --input data/sylph_phyloseq.rds --output-dir data/sylph_report
```

What it does:
- Runs `phylochart all` with your input `.rds`.
- Exports all intermediate tables into `data/sylph_report/stats/`.
- Builds the final report at `data/sylph_report/phylochart.html`.
