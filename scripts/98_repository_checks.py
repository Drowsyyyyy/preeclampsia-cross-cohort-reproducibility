"""Lightweight checks for the public portfolio checkout.

Unlike ``99_quality_checks.py``, this script does not require downloaded raw
data or the large processed expression matrices that are excluded from Git.
"""

from __future__ import annotations

import ast
import csv
import gzip
import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)
    print(f"PASS: {message}")


def read_csv_header(path: Path) -> list[str]:
    if path.suffix == ".gz":
        with gzip.open(path, mode="rt", encoding="utf-8", newline="") as handle:
            return next(csv.reader(handle))
    with path.open(mode="rt", encoding="utf-8", newline="") as handle:
        return next(csv.reader(handle))


def main() -> None:
    required = [
        ROOT / "README.md",
        ROOT / "CITATION.cff",
        ROOT / "LICENSE",
        ROOT / "config" / "analysis_plan.yaml",
        ROOT / "config" / "source_urls.yaml",
        ROOT / "data" / "processed" / "cohort_manifest.csv",
        ROOT / "data" / "processed" / "data_manifest.json",
        ROOT / "docs" / "METHODS.md",
        ROOT / "results" / "analysis_summary.json",
        ROOT / "results" / "tables" / "independent_pathway_meta_analysis.csv",
        ROOT / "results" / "tables" / "negative_control_summary.csv",
        ROOT / "output" / "pdf" / "preeclampsia_cross_cohort_report.pdf",
    ]
    for path in required:
        require(path.exists() and path.stat().st_size > 0, f"retained artifact exists: {path.relative_to(ROOT)}")

    public_text_suffixes = {".md", ".py", ".yaml", ".yml", ".cff", ".txt"}
    korean = re.compile(r"[\uac00-\ud7a3]")
    excluded = {".git", ".venv", "tmp"}
    for path in ROOT.rglob("*"):
        if not path.is_file() or path.suffix.lower() not in public_text_suffixes:
            continue
        if any(part in excluded for part in path.parts):
            continue
        source = path.read_text(encoding="utf-8")
        require(not korean.search(source), f"English-only public text: {path.relative_to(ROOT)}")
        if path.suffix == ".py":
            ast.parse(source, filename=str(path))
    print("PASS: all retained Python sources parse successfully")

    summary = json.loads((ROOT / "results" / "analysis_summary.json").read_text(encoding="utf-8"))
    require(summary["independent_cohorts"]["GSE75010_BioBank157"]["n_samples"] == 157, "discovery sample count is 157")
    require(summary["independent_meta_significant_pathways_fdr_0_05"] == 5, "retained significant-pathway count is 5")

    pathway_header = read_csv_header(ROOT / "results" / "tables" / "independent_pathway_meta_analysis.csv")
    require({"feature", "estimate", "std_error", "p_value", "fdr", "k"}.issubset(pathway_header), "pathway table exposes required statistics")

    control_header = read_csv_header(ROOT / "results" / "tables" / "negative_control_iterations.csv.gz")
    require({"level_0", "gene_set", "iteration", "estimate"}.issubset(control_header), "negative-control iterations are retained")

    for number in range(1, 7):
        matches = list((ROOT / "results" / "figures").glob(f"figure{number}_*.png"))
        require(len(matches) == 1 and matches[0].stat().st_size > 10_000, f"figure {number} is retained and nonempty")
        require(matches[0].read_bytes()[:8] == b"\x89PNG\r\n\x1a\n", f"figure {number} has a valid PNG signature")

    report = ROOT / "output" / "pdf" / "preeclampsia_cross_cohort_report.pdf"
    require(report.read_bytes()[:5] == b"%PDF-", "report has a valid PDF signature")
    print("ALL REPOSITORY CHECKS PASSED")


if __name__ == "__main__":
    main()
