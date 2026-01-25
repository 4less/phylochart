#!/usr/bin/env python3
from pathlib import Path
import argparse


TEMPLATE_PATH = Path("phylochart_template.html")
PAGES_DIR = Path("templates/pages")
DEFAULT_DATA_PATH = Path("data/taxa.tsv")
OUTPUT_PATH = Path("phylochart.html")
PLACEHOLDER = "<!-- TAXA_DATA -->"
PAGES_PLACEHOLDER = "<!-- PAGE_SECTIONS -->"
ALPHA_PLACEHOLDER = "<!-- ALPHA_DATA -->"
ALPHA_STATS_PLACEHOLDER = "<!-- ALPHA_STATS -->"
BETA_PLACEHOLDER = "<!-- BETA_DATA -->"
BETA_STATS_PLACEHOLDER = "<!-- BETA_STATS -->"
PCOA_PLACEHOLDER = "<!-- PCOA_DATA -->"

PAGE_ORDER = [
    "taxprofile",
    "alphadiv",
    "betadiv",
    "pcoa",
]


def load_page_sections() -> str:
    sections = []
    for page in PAGE_ORDER:
        page_path = PAGES_DIR / f"{page}.html"
        if not page_path.exists():
            raise SystemExit(f"Missing page template {page_path}.")
        sections.append(page_path.read_text().rstrip())
    return "\n\n".join(sections) + "\n"


def infer_alpha_path(taxa_path: Path) -> Path:
    name = taxa_path.name
    if "taxa" in name:
        return taxa_path.with_name(name.replace("taxa", "alpha", 1))
    return taxa_path.with_name("alpha.tsv")


def infer_alpha_stats_path(alpha_path: Path) -> Path:
    name = alpha_path.name
    if "alpha" in name:
        return alpha_path.with_name(name.replace("alpha", "alpha_stats", 1))
    return alpha_path.with_name("alpha_stats.tsv")


def infer_beta_path(taxa_path: Path) -> Path:
    name = taxa_path.name
    if "taxa" in name:
        return taxa_path.with_name(name.replace("taxa", "beta", 1))
    return taxa_path.with_name("beta.tsv")


def infer_beta_stats_path(beta_path: Path) -> Path:
    name = beta_path.name
    if "beta" in name:
        return beta_path.with_name(name.replace("beta", "beta_stats", 1))
    return beta_path.with_name("beta_stats.tsv")

def infer_pcoa_path(taxa_path: Path) -> Path:
    name = taxa_path.name
    if "taxa" in name:
        return taxa_path.with_name(name.replace("taxa", "pcoa", 1))
    return taxa_path.with_name("pcoa.tsv")


def split_taxa_table(text: str) -> tuple[str, str]:
    lines = text.rstrip("\n").splitlines()
    if not lines:
        return "", ""
    header = lines[0].split("\t")
    base_columns = ["Sample", "Taxon", "Rank", "Abundance", "Type"]
    base_indices = {col: header.index(col) for col in base_columns if col in header}
    if len(base_indices) != len(base_columns):
        missing = [col for col in base_columns if col not in base_indices]
        raise SystemExit(f"Missing columns in taxa TSV: {', '.join(missing)}")

    meta_columns = [col for col in header if col not in base_columns]
    sample_index = base_indices["Sample"]
    meta_rows = {}
    base_lines = ["\t".join(base_columns)]
    for raw in lines[1:]:
        if not raw.strip():
            continue
        cols = raw.split("\t")
        if len(cols) < len(header):
            cols += [""] * (len(header) - len(cols))
        base_lines.append(
            "\t".join(cols[base_indices[col]] for col in base_columns)
        )
        if meta_columns:
            sample = cols[sample_index]
            if sample and sample not in meta_rows:
                meta_rows[sample] = [cols[header.index(col)] for col in meta_columns]

    meta_lines = []
    if meta_columns:
        meta_lines.append("\t".join(["Sample", *meta_columns]))
        for sample, values in meta_rows.items():
            meta_lines.append("\t".join([sample, *values]))

    return "\n".join(base_lines), "\n".join(meta_lines)


def strip_alpha_table(text: str) -> str:
    lines = text.rstrip("\n").splitlines()
    if not lines:
        return ""
    header = lines[0].split("\t")
    base_columns = ["Sample", "Rank", "Type", "Value"]
    base_indices = {col: header.index(col) for col in base_columns if col in header}
    if len(base_indices) != len(base_columns):
        missing = [col for col in base_columns if col not in base_indices]
        raise SystemExit(f"Missing columns in alpha TSV: {', '.join(missing)}")

    base_lines = ["\t".join(base_columns)]
    for raw in lines[1:]:
        if not raw.strip():
            continue
        cols = raw.split("\t")
        if len(cols) < len(header):
            cols += [""] * (len(header) - len(cols))
        base_lines.append(
            "\t".join(cols[base_indices[col]] for col in base_columns)
        )
    return "\n".join(base_lines)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build phylochart HTML with embedded taxa/alpha data."
    )
    parser.add_argument(
        "taxa_path",
        nargs="?",
        default=None,
        help="Path to taxa TSV data (defaults to data/taxa.tsv).",
    )
    parser.add_argument(
        "alpha_path",
        nargs="?",
        default=None,
        help="Path to alpha TSV data (defaults to inferred alpha TSV).",
    )
    parser.add_argument(
        "--taxa",
        dest="taxa_flag",
        default=None,
        help="Path to taxa TSV data (overrides positional).",
    )
    parser.add_argument(
        "--alpha",
        dest="alpha_flag",
        default=None,
        help="Path to alpha TSV data (overrides positional).",
    )
    parser.add_argument(
        "--alpha-stats",
        dest="alpha_stats_flag",
        default=None,
        help="Path to alpha stats TSV data (overrides inferred path).",
    )
    parser.add_argument(
        "--beta",
        dest="beta_flag",
        default=None,
        help="Path to beta diversity TSV data.",
    )
    parser.add_argument(
        "--pcoa",
        dest="pcoa_flag",
        default=None,
        help="Placeholder path to PCoA data (not yet used).",
    )
    parser.add_argument(
        "--beta-stats",
        dest="beta_stats_flag",
        default=None,
        help="Path to beta stats TSV data (overrides inferred path).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.pcoa_flag:
        # Placeholder flag for future PCoA embedding support.
        pass

    data_path = (
        Path(args.taxa_flag)
        if args.taxa_flag
        else Path(args.taxa_path)
        if args.taxa_path
        else DEFAULT_DATA_PATH
    )
    alpha_path = (
        Path(args.alpha_flag)
        if args.alpha_flag
        else Path(args.alpha_path)
        if args.alpha_path
        else None
    )
    alpha_stats_path = Path(args.alpha_stats_flag) if args.alpha_stats_flag else None
    beta_path = Path(args.beta_flag) if args.beta_flag else None
    beta_stats_path = Path(args.beta_stats_flag) if args.beta_stats_flag else None
    template = TEMPLATE_PATH.read_text()
    if PLACEHOLDER not in template:
        raise SystemExit(f"Missing placeholder {PLACEHOLDER} in {TEMPLATE_PATH}.")
    if PAGES_PLACEHOLDER not in template:
        raise SystemExit(f"Missing placeholder {PAGES_PLACEHOLDER} in {TEMPLATE_PATH}.")
    if ALPHA_PLACEHOLDER not in template:
        raise SystemExit(f"Missing placeholder {ALPHA_PLACEHOLDER} in {TEMPLATE_PATH}.")
    if ALPHA_STATS_PLACEHOLDER not in template:
        raise SystemExit(
            f"Missing placeholder {ALPHA_STATS_PLACEHOLDER} in {TEMPLATE_PATH}."
        )
    if BETA_PLACEHOLDER not in template:
        raise SystemExit(f"Missing placeholder {BETA_PLACEHOLDER} in {TEMPLATE_PATH}.")
    if BETA_STATS_PLACEHOLDER not in template:
        raise SystemExit(
            f"Missing placeholder {BETA_STATS_PLACEHOLDER} in {TEMPLATE_PATH}."
        )
    if PCOA_PLACEHOLDER not in template:
        raise SystemExit(
            f"Missing placeholder {PCOA_PLACEHOLDER} in {TEMPLATE_PATH}."
        )

    pages = load_page_sections()

    data = data_path.read_text().rstrip("\n")
    taxa_base, taxa_meta = split_taxa_table(data)
    embedded = f'\n    <script type="text/plain" id="taxaData">\n{taxa_base}\n    </script>\n'
    if taxa_meta:
        embedded += f'\n    <script type="text/plain" id="taxaMeta">\n{taxa_meta}\n    </script>\n'

    if alpha_path is None:
        inferred = infer_alpha_path(data_path)
        if inferred.exists():
            alpha_path = inferred
    alpha_embedded = "\n"
    if alpha_path and alpha_path.exists():
        alpha_text = alpha_path.read_text().rstrip("\n")
        alpha_base = strip_alpha_table(alpha_text)
        alpha_embedded = (
            f'\n    <script type="text/plain" id="alphaData">\n{alpha_base}\n    </script>\n'
        )

    alpha_stats_embedded = "\n"
    if alpha_stats_path is None and alpha_path:
        inferred_stats = infer_alpha_stats_path(alpha_path)
        if inferred_stats.exists():
            alpha_stats_path = inferred_stats
    if alpha_stats_path and alpha_stats_path.exists():
        stats_text = alpha_stats_path.read_text().rstrip("\n")
        alpha_stats_embedded = (
            f'\n    <script type="text/plain" id="alphaStats">\n{stats_text}\n    </script>\n'
        )

    if beta_path is None:
        inferred = infer_beta_path(data_path)
        if inferred.exists():
            beta_path = inferred
    beta_embedded = "\n"
    if beta_path and beta_path.exists():
        beta_text = beta_path.read_text().rstrip("\n")
        beta_embedded = (
            f'\n    <script type="text/plain" id="betaData">\n{beta_text}\n    </script>\n'
        )

    beta_stats_embedded = "\n"
    if beta_stats_path is None and beta_path:
        inferred_stats = infer_beta_stats_path(beta_path)
        if inferred_stats.exists():
            beta_stats_path = inferred_stats
    if beta_stats_path and beta_stats_path.exists():
        stats_text = beta_stats_path.read_text().rstrip("\n")
        beta_stats_embedded = (
            f'\n    <script type="text/plain" id="betaStats">\n{stats_text}\n    </script>\n'
        )

    pcoa_path = Path(args.pcoa_flag) if args.pcoa_flag else None
    if pcoa_path is None:
        inferred = infer_pcoa_path(data_path)
        if inferred.exists():
            pcoa_path = inferred
    pcoa_embedded = "\n"
    if pcoa_path and pcoa_path.exists():
        pcoa_text = pcoa_path.read_text().rstrip("\n")
        pcoa_embedded = (
            f'\n    <script type="text/plain" id="pcoaData">\n{pcoa_text}\n    </script>\n'
        )

    output = template.replace(PAGES_PLACEHOLDER, pages, 1)
    output = output.replace(PLACEHOLDER, embedded, 1)
    output = output.replace(ALPHA_PLACEHOLDER, alpha_embedded, 1)
    output = output.replace(ALPHA_STATS_PLACEHOLDER, alpha_stats_embedded, 1)
    output = output.replace(BETA_PLACEHOLDER, beta_embedded, 1)
    output = output.replace(BETA_STATS_PLACEHOLDER, beta_stats_embedded, 1)
    output = output.replace(PCOA_PLACEHOLDER, pcoa_embedded, 1)
    OUTPUT_PATH.write_text(output)
    print(f"Wrote {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
