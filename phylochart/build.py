from __future__ import annotations

import argparse
import shutil
from pathlib import Path

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


def _package_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _resolve_template_path(template_path: Path | None) -> Path:
    if template_path is not None:
        return template_path
    for candidate in (
        Path.cwd() / "phylochart_template.html",
        _package_root() / "phylochart_template.html",
    ):
        if candidate.exists():
            return candidate
    raise SystemExit("Could not find phylochart_template.html.")


def _resolve_pages_dir(pages_dir: Path | None) -> Path:
    if pages_dir is not None:
        return pages_dir
    for candidate in (
        Path.cwd() / "templates/pages",
        _package_root() / "templates/pages",
    ):
        if candidate.exists():
            return candidate
    raise SystemExit("Could not find templates/pages directory.")


def _resolve_assets_dir() -> Path | None:
    for candidate in (
        Path.cwd() / "assets",
        _package_root() / "assets",
    ):
        if candidate.exists():
            return candidate
    return None


def _copy_assets_for_report(output_path: Path) -> None:
    assets_src = _resolve_assets_dir()
    if assets_src is None:
        return
    output_dir = output_path.parent
    output_dir.mkdir(parents=True, exist_ok=True)
    assets_dst = output_dir / "assets"
    if assets_src.resolve() == assets_dst.resolve():
        return
    shutil.copytree(assets_src, assets_dst, dirs_exist_ok=True)


def load_page_sections(pages_dir: Path) -> str:
    sections = []
    for page in PAGE_ORDER:
        page_path = pages_dir / f"{page}.html"
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
    meta_rows: dict[str, list[str]] = {}
    base_lines = ["\t".join(base_columns)]
    for raw in lines[1:]:
        if not raw.strip():
            continue
        cols = raw.split("\t")
        if len(cols) < len(header):
            cols += [""] * (len(header) - len(cols))
        base_lines.append("\t".join(cols[base_indices[col]] for col in base_columns))
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
        base_lines.append("\t".join(cols[base_indices[col]] for col in base_columns))
    return "\n".join(base_lines)


def build_html(
    taxa_path: Path,
    output_path: Path,
    *,
    alpha_path: Path | None = None,
    alpha_stats_path: Path | None = None,
    beta_path: Path | None = None,
    beta_stats_path: Path | None = None,
    pcoa_path: Path | None = None,
    template_path: Path | None = None,
    pages_dir: Path | None = None,
) -> Path:
    template_path = _resolve_template_path(template_path)
    pages_dir = _resolve_pages_dir(pages_dir)
    template = template_path.read_text()
    for placeholder in (
        PLACEHOLDER,
        PAGES_PLACEHOLDER,
        ALPHA_PLACEHOLDER,
        ALPHA_STATS_PLACEHOLDER,
        BETA_PLACEHOLDER,
        BETA_STATS_PLACEHOLDER,
        PCOA_PLACEHOLDER,
    ):
        if placeholder not in template:
            raise SystemExit(f"Missing placeholder {placeholder} in {template_path}.")

    pages = load_page_sections(pages_dir)

    data = taxa_path.read_text().rstrip("\n")
    taxa_base, taxa_meta = split_taxa_table(data)
    taxa_embedded = f'\n    <script type="text/plain" id="taxaData">\n{taxa_base}\n    </script>\n'
    if taxa_meta:
        taxa_embedded += (
            f'\n    <script type="text/plain" id="taxaMeta">\n{taxa_meta}\n    </script>\n'
        )

    if alpha_path is None:
        inferred = infer_alpha_path(taxa_path)
        if inferred.exists():
            alpha_path = inferred
    alpha_embedded = "\n"
    if alpha_path and alpha_path.exists():
        alpha_base = strip_alpha_table(alpha_path.read_text().rstrip("\n"))
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
        inferred = infer_beta_path(taxa_path)
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

    if pcoa_path is None:
        inferred = infer_pcoa_path(taxa_path)
        if inferred.exists():
            pcoa_path = inferred
    pcoa_embedded = "\n"
    if pcoa_path and pcoa_path.exists():
        pcoa_text = pcoa_path.read_text().rstrip("\n")
        pcoa_embedded = (
            f'\n    <script type="text/plain" id="pcoaData">\n{pcoa_text}\n    </script>\n'
        )

    output = template.replace(PAGES_PLACEHOLDER, pages, 1)
    output = output.replace(PLACEHOLDER, taxa_embedded, 1)
    output = output.replace(ALPHA_PLACEHOLDER, alpha_embedded, 1)
    output = output.replace(ALPHA_STATS_PLACEHOLDER, alpha_stats_embedded, 1)
    output = output.replace(BETA_PLACEHOLDER, beta_embedded, 1)
    output = output.replace(BETA_STATS_PLACEHOLDER, beta_stats_embedded, 1)
    output = output.replace(PCOA_PLACEHOLDER, pcoa_embedded, 1)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(output)
    _copy_assets_for_report(output_path)
    return output_path


def _parse_build_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build phylochart HTML with embedded taxa/alpha/beta/pcoa data.",
    )
    parser.add_argument(
        "taxa_path",
        nargs="?",
        default=None,
        help="Path to taxa TSV data (defaults to data/taxa.tsv).",
    )
    parser.add_argument("--taxa", dest="taxa_flag", default=None, help="Path to taxa TSV.")
    parser.add_argument("--alpha", dest="alpha_flag", default=None, help="Path to alpha TSV.")
    parser.add_argument(
        "--alpha-stats",
        dest="alpha_stats_flag",
        default=None,
        help="Path to alpha stats TSV.",
    )
    parser.add_argument("--beta", dest="beta_flag", default=None, help="Path to beta TSV.")
    parser.add_argument(
        "--beta-stats",
        dest="beta_stats_flag",
        default=None,
        help="Path to beta stats TSV.",
    )
    parser.add_argument("--pcoa", dest="pcoa_flag", default=None, help="Path to pcoa TSV.")
    parser.add_argument(
        "--output",
        default="phylochart.html",
        help="Output HTML path (default: phylochart.html).",
    )
    parser.add_argument(
        "--template",
        default=None,
        help="Optional template HTML path (default: phylochart_template.html).",
    )
    parser.add_argument(
        "--pages-dir",
        default=None,
        help="Optional pages dir path (default: templates/pages).",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = _parse_build_args(argv)
    taxa_path = Path(args.taxa_flag) if args.taxa_flag else Path(args.taxa_path or "data/taxa.tsv")
    output_path = Path(args.output)
    path = build_html(
        taxa_path=taxa_path,
        output_path=output_path,
        alpha_path=Path(args.alpha_flag) if args.alpha_flag else None,
        alpha_stats_path=Path(args.alpha_stats_flag) if args.alpha_stats_flag else None,
        beta_path=Path(args.beta_flag) if args.beta_flag else None,
        beta_stats_path=Path(args.beta_stats_flag) if args.beta_stats_flag else None,
        pcoa_path=Path(args.pcoa_flag) if args.pcoa_flag else None,
        template_path=Path(args.template) if args.template else None,
        pages_dir=Path(args.pages_dir) if args.pages_dir else None,
    )
    print(f"Wrote {path}")
    return 0
