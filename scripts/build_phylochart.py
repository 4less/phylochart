#!/usr/bin/env python3
from pathlib import Path
import sys


TEMPLATE_PATH = Path("phylochart_template.html")
DEFAULT_DATA_PATH = Path("data/taxa.tsv")
OUTPUT_PATH = Path("phylochart.html")
PLACEHOLDER = "<!-- TAXA_DATA -->"


def main() -> None:
    data_path = DEFAULT_DATA_PATH
    if len(sys.argv) > 1:
        data_path = Path(sys.argv[1])
    template = TEMPLATE_PATH.read_text()
    if PLACEHOLDER not in template:
        raise SystemExit(f"Missing placeholder {PLACEHOLDER} in {TEMPLATE_PATH}.")

    data = data_path.read_text().rstrip("\n")
    embedded = f'\n    <script type="text/plain" id="taxaData">\n{data}\n    </script>\n'
    output = template.replace(PLACEHOLDER, embedded, 1)
    OUTPUT_PATH.write_text(output)
    print(f"Wrote {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
