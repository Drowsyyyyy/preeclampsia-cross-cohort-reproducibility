"""Fail-fast integrity checks for retained data and analysis artifacts."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
PROCESSED = ROOT / "data" / "processed"
TABLES = ROOT / "results" / "tables"
FIGURES = ROOT / "results" / "figures"


def file_hash(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)
    print(f"PASS: {message}")


def main() -> None:
    manifest = json.loads((PROCESSED / "data_manifest.json").read_text(encoding="utf-8"))
    for item in manifest["files"]:
        path = ROOT / item["path"]
        require(path.exists(), f"source file exists: {item['path']}")
        require(path.stat().st_size == item["bytes"], f"source byte count matches: {item['path']}")
        require(file_hash(path) == item["sha256"], f"source SHA-256 matches: {item['path']}")

    cohorts = pd.read_csv(PROCESSED / "cohort_manifest.csv")
    primary = cohorts.loc[cohorts["included_primary"].astype(bool)]

    expected = {
        ("GSE75010_BioBank157", "PE"): 80,
        ("GSE75010_BioBank157", "Control"): 77,
        ("GSE190971", "PE"): 7,
        ("GSE190971", "Control"): 6,
        ("GSE204835", "PE"): 12,
        ("GSE204835", "Control"): 12,
    }
    counts = primary.groupby(["cohort_id", "group"]).size().to_dict()
    for key, value in expected.items():
        require(counts.get(key) == value, f"primary sample count {key[0]} {key[1]} = {value}")

    historical_ids = {
        "GSE30186",
        "GSE10588",
        "GSE24129",
        "GSE25906",
        "GSE43942",
        "GSE4707",
        "GSE44711",
    }
    historical_n = primary.loc[primary["cohort_id"].isin(historical_ids)].shape[0]
    require(historical_n == 173, "historical embedded PE/control samples total 173")
    require(
        not primary.loc[primary["cohort_id"].isin(historical_ids), "independently_collected"].astype(bool).any(),
        "historical embedded samples are never marked independently collected",
    )

    expression_files = {
        "GSE75010_BioBank157": "GSE75010_BioBank157_expression.csv.gz",
        "GSE75010_historical173": "GSE75010_historical173_expression.csv.gz",
        "GSE190971": "GSE190971_expression.csv.gz",
        "GSE204835": "GSE204835_expression.csv.gz",
        "GSE177049_gene": "GSE177049_gene_expression.csv.gz",
        "GSE177049_mirna": "GSE177049_mirna_expression.csv.gz",
    }
    matrices: dict[str, pd.DataFrame] = {}
    for key, name in expression_files.items():
        matrix = pd.read_csv(PROCESSED / name, index_col=0)
        matrices[key] = matrix
        require(matrix.index.is_unique, f"{key} feature identifiers are unique")
        require(matrix.columns.is_unique, f"{key} sample identifiers are unique")
        require(np.isfinite(matrix.to_numpy(dtype=float)).all(), f"{key} matrix contains only finite values")

    require(matrices["GSE75010_BioBank157"].shape[1] == 157, "discovery matrix has 157 samples")
    require(matrices["GSE75010_historical173"].shape[1] == 173, "historical matrix has 173 samples")
    require(
        set(matrices["GSE75010_BioBank157"].columns).isdisjoint(matrices["GSE75010_historical173"].columns),
        "discovery and historical matrix columns are disjoint",
    )

    meta = pd.read_csv(TABLES / "independent_pathway_meta_analysis.csv")
    require(meta["k"].between(1, 3).all(), "pathway meta-analysis uses at most three independent cohorts")
    require("cohort" not in meta.columns, "pathway meta table contains pooled effects rather than raw-scale pooled samples")

    target = pd.read_csv(TABLES / "target_gene_cohort_effects.csv")
    require(
        set(target["cohort"]).issubset({"GSE75010_BioBank157", "GSE190971", "GSE204835_term"}),
        "independent target-gene analysis excludes embedded historical studies",
    )

    controls = pd.read_csv(TABLES / "negative_control_iterations.csv.gz")
    control_counts = controls.groupby(["gene_set", "level_0"]).size()
    require((control_counts.xs("matched_random", level="level_0") == 500).all(), "500 matched-random controls per measured gene set")
    require((control_counts.xs("label_permutation", level="level_0") == 500).all(), "500 label permutations per measured gene set")

    required_tables = [
        "discovery_gene_differential_expression.csv.gz",
        "independent_gene_meta_analysis.csv.gz",
        "independent_pathway_meta_analysis.csv",
        "leave_one_cohort_out_meta_analysis.csv.gz",
        "composition_adjusted_pathways.csv",
        "negative_control_summary.csv",
        "GSE177049_mirna_effects.csv",
    ]
    for name in required_tables:
        require((TABLES / name).exists() and (TABLES / name).stat().st_size > 0, f"analysis table retained: {name}")
    for number in range(1, 7):
        matches = list(FIGURES.glob(f"figure{number}_*.png"))
        require(len(matches) == 1 and matches[0].stat().st_size > 10_000, f"figure {number} retained and nonempty")

    summary = json.loads((ROOT / "results" / "analysis_summary.json").read_text(encoding="utf-8"))
    require(summary["independent_cohorts"]["GSE75010_BioBank157"]["n_samples"] == 157, "summary discovery n is internally consistent")
    require(summary["independent_meta_significant_pathways_fdr_0_05"] == 5, "summary pathway count matches retained analysis")
    print("ALL QUALITY CHECKS PASSED")


if __name__ == "__main__":
    main()
