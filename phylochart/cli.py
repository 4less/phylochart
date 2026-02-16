from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

from .build import build_html

R_SCRIPT_BY_COMMAND = {
    "taxa": "export_ps_taxa.R",
    "alpha": "export_ps_alpha.R",
    "alpha_stats": "export_ps_alpha_stats.R",
    "beta": "export_ps_beta.R",
    "beta_stats": "export_ps_beta_stats.R",
    "pcoa": "export_ps_pcoa.R",
}


def _package_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _resolve_r_script(script_name: str) -> Path:
    for candidate in (
        Path.cwd() / "src" / script_name,
        _package_root() / "src" / script_name,
    ):
        if candidate.exists():
            return candidate
    raise SystemExit(f"Could not find R script: {script_name}")


def _run_r_script(script_name: str, args: list[str]) -> None:
    script_path = _resolve_r_script(script_name)
    cmd = ["Rscript", str(script_path), *args]
    subprocess.run(cmd, check=True)


def _infer_stats_output(output_path: Path, token: str, default_name: str) -> Path:
    name = output_path.name
    if token in name:
        return output_path.with_name(name.replace(token, f"{token}_stats", 1))
    return output_path.with_name(default_name)


def _filter_args(filter_treatment: bool) -> list[str]:
    if filter_treatment:
        return ["--filter-treatment"]
    return ["--no-filter-treatment"]


def _add_filter_flags(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--filter-treatment",
        dest="filter_treatment",
        action="store_true",
        default=True,
        help="Filter rows with missing Treatment metadata (default).",
    )
    parser.add_argument(
        "--no-filter-treatment",
        dest="filter_treatment",
        action="store_false",
        help="Do not filter missing Treatment metadata.",
    )


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="phylochart command line interface")
    subparsers = parser.add_subparsers(dest="command", required=True)

    taxa = subparsers.add_parser("taxa", help="Export taxa table from a phyloseq RDS.")
    taxa.add_argument("--input", required=True, help="Input phyloseq RDS path.")
    taxa.add_argument("--output", required=True, help="Output taxa TSV path.")
    taxa.add_argument(
        "--types",
        choices=("relative", "absolute", "both"),
        default="both",
        help="Abundance types to export (default: both).",
    )
    taxa.add_argument(
        "--drop-zeros",
        dest="drop_zeros",
        action="store_true",
        default=True,
        help="Drop zero-abundance rows (default).",
    )
    taxa.add_argument(
        "--keep-zeros",
        dest="drop_zeros",
        action="store_false",
        help="Keep zero-abundance rows.",
    )
    _add_filter_flags(taxa)

    alpha = subparsers.add_parser("alpha", help="Export alpha diversity and stats.")
    alpha.add_argument("--input", required=True, help="Input phyloseq RDS path.")
    alpha.add_argument("--output", required=True, help="Output alpha TSV path.")
    alpha.add_argument(
        "--output-stats",
        default=None,
        help="Output alpha stats TSV path (default: inferred from --output).",
    )
    _add_filter_flags(alpha)

    alpha_stats = subparsers.add_parser("alpha-stats", help="Export alpha diversity stats only.")
    alpha_stats.add_argument("--input", required=True, help="Input phyloseq RDS path.")
    alpha_stats.add_argument("--output", required=True, help="Output alpha stats TSV path.")
    _add_filter_flags(alpha_stats)

    beta = subparsers.add_parser("beta", help="Export beta diversity and stats.")
    beta.add_argument("--input", required=True, help="Input phyloseq RDS path.")
    beta.add_argument("--output", required=True, help="Output beta TSV path.")
    beta.add_argument(
        "--output-stats",
        default=None,
        help="Output beta stats TSV path (default: inferred from --output).",
    )
    _add_filter_flags(beta)

    beta_stats = subparsers.add_parser("beta-stats", help="Export beta diversity stats only.")
    beta_stats.add_argument("--input", required=True, help="Input phyloseq RDS path.")
    beta_stats.add_argument("--output", required=True, help="Output beta stats TSV path.")
    _add_filter_flags(beta_stats)

    pcoa = subparsers.add_parser("pcoa", help="Export PCoA/NMDS coordinates.")
    pcoa.add_argument("--input", required=True, help="Input phyloseq RDS path.")
    pcoa.add_argument("--output", required=True, help="Output PCoA TSV path.")
    _add_filter_flags(pcoa)

    all_cmd = subparsers.add_parser(
        "all",
        help="Run taxa/alpha/beta/pcoa exports and build phylochart report.",
    )
    all_cmd.add_argument("--input", required=True, help="Input phyloseq RDS path.")
    all_cmd.add_argument(
        "--output-dir",
        required=True,
        help="Output directory. Writes stats under output-dir/stats and report at output-dir/phylochart.html.",
    )
    all_cmd.add_argument(
        "--types",
        choices=("relative", "absolute", "both"),
        default="both",
        help="Abundance types for taxa export (default: both).",
    )
    all_cmd.add_argument(
        "--drop-zeros",
        dest="drop_zeros",
        action="store_true",
        default=True,
        help="Drop zero-abundance rows for taxa export (default).",
    )
    all_cmd.add_argument(
        "--keep-zeros",
        dest="drop_zeros",
        action="store_false",
        help="Keep zero-abundance rows for taxa export.",
    )
    all_cmd.add_argument("--template", default=None, help="Optional template path.")
    all_cmd.add_argument("--pages-dir", default=None, help="Optional pages directory path.")
    _add_filter_flags(all_cmd)

    build = subparsers.add_parser("build", help="Build final phylochart HTML.")
    build.add_argument("--taxa", required=True, help="Input taxa TSV path.")
    build.add_argument("--alpha", default=None, help="Input alpha TSV path.")
    build.add_argument("--alpha-stats", default=None, help="Input alpha stats TSV path.")
    build.add_argument("--beta", default=None, help="Input beta TSV path.")
    build.add_argument("--beta-stats", default=None, help="Input beta stats TSV path.")
    build.add_argument("--pcoa", default=None, help="Input pcoa TSV path.")
    build.add_argument(
        "--output",
        default="phylochart.html",
        help="Output HTML path (default: phylochart.html).",
    )
    build.add_argument("--template", default=None, help="Optional template path.")
    build.add_argument("--pages-dir", default=None, help="Optional pages directory path.")

    return parser


def _run_taxa(args: argparse.Namespace) -> int:
    cmd = [
        "--input",
        args.input,
        "--output",
        args.output,
        "--types",
        args.types,
        *(_filter_args(args.filter_treatment)),
        "--drop-zeros" if args.drop_zeros else "--keep-zeros",
    ]
    _run_r_script(R_SCRIPT_BY_COMMAND["taxa"], cmd)
    return 0


def _run_alpha(args: argparse.Namespace) -> int:
    output = Path(args.output)
    output_stats = Path(args.output_stats) if args.output_stats else _infer_stats_output(
        output, "alpha", "alpha_stats.tsv"
    )
    _run_r_script(
        R_SCRIPT_BY_COMMAND["alpha"],
        ["--input", args.input, "--output", str(output), *(_filter_args(args.filter_treatment))],
    )
    _run_r_script(
        R_SCRIPT_BY_COMMAND["alpha_stats"],
        [
            "--input",
            args.input,
            "--output",
            str(output_stats),
            *(_filter_args(args.filter_treatment)),
        ],
    )
    return 0


def _run_beta(args: argparse.Namespace) -> int:
    output = Path(args.output)
    output_stats = Path(args.output_stats) if args.output_stats else _infer_stats_output(
        output, "beta", "beta_stats.tsv"
    )
    _run_r_script(
        R_SCRIPT_BY_COMMAND["beta"],
        ["--input", args.input, "--output", str(output), *(_filter_args(args.filter_treatment))],
    )
    _run_r_script(
        R_SCRIPT_BY_COMMAND["beta_stats"],
        [
            "--input",
            args.input,
            "--output",
            str(output_stats),
            *(_filter_args(args.filter_treatment)),
        ],
    )
    return 0


def _run_alpha_stats(args: argparse.Namespace) -> int:
    _run_r_script(
        R_SCRIPT_BY_COMMAND["alpha_stats"],
        ["--input", args.input, "--output", args.output, *(_filter_args(args.filter_treatment))],
    )
    return 0


def _run_beta_stats(args: argparse.Namespace) -> int:
    _run_r_script(
        R_SCRIPT_BY_COMMAND["beta_stats"],
        ["--input", args.input, "--output", args.output, *(_filter_args(args.filter_treatment))],
    )
    return 0


def _run_pcoa(args: argparse.Namespace) -> int:
    _run_r_script(
        R_SCRIPT_BY_COMMAND["pcoa"],
        ["--input", args.input, "--output", args.output, *(_filter_args(args.filter_treatment))],
    )
    return 0


def _run_build(args: argparse.Namespace) -> int:
    output_path = Path(args.output)
    path = build_html(
        taxa_path=Path(args.taxa),
        output_path=output_path,
        alpha_path=Path(args.alpha) if args.alpha else None,
        alpha_stats_path=Path(args.alpha_stats) if args.alpha_stats else None,
        beta_path=Path(args.beta) if args.beta else None,
        beta_stats_path=Path(args.beta_stats) if args.beta_stats else None,
        pcoa_path=Path(args.pcoa) if args.pcoa else None,
        template_path=Path(args.template) if args.template else None,
        pages_dir=Path(args.pages_dir) if args.pages_dir else None,
    )
    print(f"Wrote {path}")
    return 0


def _run_all(args: argparse.Namespace) -> int:
    output_dir = Path(args.output_dir)
    stats_dir = output_dir / "stats"
    stats_dir.mkdir(parents=True, exist_ok=True)

    taxa_path = stats_dir / "taxa.tsv"
    alpha_path = stats_dir / "alpha.tsv"
    alpha_stats_path = stats_dir / "alpha_stats.tsv"
    beta_path = stats_dir / "beta.tsv"
    beta_stats_path = stats_dir / "beta_stats.tsv"
    pcoa_path = stats_dir / "pcoa.tsv"
    report_path = output_dir / "phylochart.html"

    _run_r_script(
        R_SCRIPT_BY_COMMAND["taxa"],
        [
            "--input",
            args.input,
            "--output",
            str(taxa_path),
            "--types",
            args.types,
            *(_filter_args(args.filter_treatment)),
            "--drop-zeros" if args.drop_zeros else "--keep-zeros",
        ],
    )
    _run_r_script(
        R_SCRIPT_BY_COMMAND["alpha"],
        [
            "--input",
            args.input,
            "--output",
            str(alpha_path),
            *(_filter_args(args.filter_treatment)),
        ],
    )
    _run_r_script(
        R_SCRIPT_BY_COMMAND["alpha_stats"],
        [
            "--input",
            args.input,
            "--output",
            str(alpha_stats_path),
            *(_filter_args(args.filter_treatment)),
        ],
    )
    _run_r_script(
        R_SCRIPT_BY_COMMAND["beta"],
        [
            "--input",
            args.input,
            "--output",
            str(beta_path),
            *(_filter_args(args.filter_treatment)),
        ],
    )
    _run_r_script(
        R_SCRIPT_BY_COMMAND["beta_stats"],
        [
            "--input",
            args.input,
            "--output",
            str(beta_stats_path),
            *(_filter_args(args.filter_treatment)),
        ],
    )
    _run_r_script(
        R_SCRIPT_BY_COMMAND["pcoa"],
        [
            "--input",
            args.input,
            "--output",
            str(pcoa_path),
            *(_filter_args(args.filter_treatment)),
        ],
    )

    path = build_html(
        taxa_path=taxa_path,
        output_path=report_path,
        alpha_path=alpha_path,
        alpha_stats_path=alpha_stats_path,
        beta_path=beta_path,
        beta_stats_path=beta_stats_path,
        pcoa_path=pcoa_path,
        template_path=Path(args.template) if args.template else None,
        pages_dir=Path(args.pages_dir) if args.pages_dir else None,
    )
    print(f"Wrote {path}")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    try:
        if args.command == "taxa":
            return _run_taxa(args)
        if args.command == "alpha":
            return _run_alpha(args)
        if args.command == "alpha-stats":
            return _run_alpha_stats(args)
        if args.command == "beta":
            return _run_beta(args)
        if args.command == "beta-stats":
            return _run_beta_stats(args)
        if args.command == "pcoa":
            return _run_pcoa(args)
        if args.command == "all":
            return _run_all(args)
        if args.command == "build":
            return _run_build(args)
        parser.error(f"Unknown command: {args.command}")
    except subprocess.CalledProcessError as exc:
        return exc.returncode
    return 1


if __name__ == "__main__":
    sys.exit(main())
