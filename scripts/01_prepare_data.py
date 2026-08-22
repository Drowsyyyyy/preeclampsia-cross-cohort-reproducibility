"""Build the sample-provenance manifest and analysis-ready matrices.

This script is deliberately strict about the GSE75010 nesting problem:
the 157 newly collected RCWIH BioBank samples form the discovery cohort,
whereas samples imported from seven earlier GEO studies are retained only as
historical/internal replication cohorts. GSE190971 and GSE204835 are the two
independently collected external RNA-seq resources.
"""

from __future__ import annotations

import gzip
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

from common import (
    filter_count_matrix,
    first_description_token,
    median_ratio_log_expression,
    parse_gestational_weeks,
    parse_number,
    parse_soft_samples,
    sha256_file,
    write_json,
)


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
PROCESSED = ROOT / "data" / "processed"
PROCESSED.mkdir(parents=True, exist_ok=True)

HISTORICAL = [
    "GSE30186",
    "GSE10588",
    "GSE24129",
    "GSE25906",
    "GSE43942",
    "GSE4707",
    "GSE44711",
]


def clean_sex(value: object) -> str:
    text = str(value).strip().upper() if value is not None else ""
    if text.startswith("F"):
        return "F"
    if text.startswith("M"):
        return "M"
    return "Unknown"


def classify_historical(gse: str, sample: dict[str, object]) -> str:
    title = str(sample.get("title", "")).lower()
    source = str(sample.get("source_name_ch1", "")).lower()
    state = str(
        sample.get(
            "disease_state",
            sample.get("classification", sample.get("phenotype", sample.get("condition", ""))),
        )
    ).lower()
    characteristics = str(sample.get("characteristics_raw", "")).lower()
    combined = " | ".join([title, source, state, characteristics])
    if gse == "GSE24129" and "growth restriction" in combined:
        return "FGR_other"
    if any(token in combined for token in ["preeclampsia", "pre-eclampsia", "preeclamptic", "eopet", "pe(eo)", "pe(lo)"]):
        return "PE"
    if any(token in combined for token in ["normal", "control"]):
        return "Control"
    return "Unclear"


def prepare_gse75010() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    complete = pd.read_csv(RAW / "GSE75010_complete_dataset.csv.gz", index_col=0)
    complete.index = complete.index.astype(str).str.upper()
    complete = complete.groupby(complete.index, sort=False).mean()

    new_samples = parse_soft_samples(RAW / "GSE75010_family.soft.gz")
    new_rows: list[dict[str, object]] = []
    for sample in new_samples:
        matrix_column = first_description_token(sample.get("description"))
        ga_week = parse_number(sample.get("ga_week"))
        ga_day = parse_number(sample.get("ga_day"))
        ga = ga_week + ga_day / 7.0 if np.isfinite(ga_week) else np.nan
        diagnosis = str(sample.get("diagnosis", "")).strip().upper()
        group = "PE" if diagnosis == "PE" else "Control" if diagnosis == "NON-PE" else "Unclear"
        title = str(sample.get("title", ""))
        new_rows.append(
            {
                "cohort_id": "GSE75010_BioBank157",
                "original_study": "GSE75010_new_BioBank",
                "geo_accession": sample["geo_accession"],
                "matrix_column": matrix_column,
                "group": group,
                "pe": 1 if group == "PE" else 0 if group == "Control" else np.nan,
                "gestational_age_weeks": ga,
                "gestational_age_category": "Preterm" if np.isfinite(ga) and ga < 34 else "Term",
                "fetal_sex": clean_sex(sample.get("infant_gender")),
                "chronic_hypertension": int("-CH" in title.upper()),
                "tissue": "placenta",
                "platform": "Affymetrix Human Gene 1.0 ST Array (GPL6244)",
                "role": "discovery",
                "independently_collected": True,
                "included_primary": group in {"PE", "Control"} and matrix_column in complete.columns,
                "exclusion_reason": "" if matrix_column in complete.columns else "missing_from_integrated_matrix",
                "source_file": "GSE75010_complete_dataset.csv.gz",
                "provenance_notes": "New RCWIH BioBank sample registered directly under GSE75010",
            }
        )

    historical_rows: list[dict[str, object]] = []
    for gse in HISTORICAL:
        samples = parse_soft_samples(RAW / f"{gse}_family.soft.gz")
        for sample in samples:
            group = classify_historical(gse, sample)
            ga_value = sample.get("gestational_age", sample.get("gestational_age_weeks"))
            ga = parse_gestational_weeks(ga_value)
            sex = clean_sex(sample.get("gender", sample.get("infant_gender")))
            gsm = str(sample["geo_accession"])
            included = gsm in complete.columns and group in {"PE", "Control"}
            reason = ""
            if gsm not in complete.columns:
                reason = "not_selected_in_GSE75010_integrated_matrix"
            elif group not in {"PE", "Control"}:
                reason = "non_PE_comparator"
            historical_rows.append(
                {
                    "cohort_id": gse,
                    "original_study": gse,
                    "geo_accession": gsm,
                    "matrix_column": gsm,
                    "group": group,
                    "pe": 1 if group == "PE" else 0 if group == "Control" else np.nan,
                    "gestational_age_weeks": ga,
                    "gestational_age_category": "Preterm" if np.isfinite(ga) and ga < 37 else "Term" if np.isfinite(ga) else "Unknown",
                    "fetal_sex": sex,
                    "chronic_hypertension": np.nan,
                    "tissue": str(sample.get("source_name_ch1", "placenta")),
                    "platform": "Historical platform harmonized by GSE75010",
                    "role": "historical_internal_replication",
                    "independently_collected": False,
                    "included_primary": included,
                    "exclusion_reason": reason,
                    "source_file": "GSE75010_complete_dataset.csv.gz",
                    "provenance_notes": "Imported into GSE75010 from the named original GEO study; not external validation",
                }
            )

    new_manifest = pd.DataFrame(new_rows)
    historical_manifest = pd.DataFrame(historical_rows)
    discovery_cols = new_manifest.loc[new_manifest["included_primary"], "matrix_column"].tolist()
    historical_cols = historical_manifest.loc[historical_manifest["included_primary"], "matrix_column"].tolist()
    discovery = complete.loc[:, discovery_cols]
    historical = complete.loc[:, historical_cols]
    return new_manifest, historical_manifest, pd.concat(
        {"discovery": discovery, "historical": historical}, axis=1
    )


def prepare_gse190971() -> tuple[pd.DataFrame, pd.DataFrame]:
    raw = pd.read_csv(RAW / "GSE190971_Raw_gene_counts_matrix_PLAC.txt.gz", sep="\t")
    raw["Gene_Symbol"] = raw["Gene_Symbol"].astype(str).str.upper()
    counts = raw.groupby("Gene_Symbol", sort=False).sum(numeric_only=True)
    counts = filter_count_matrix(counts, min_count=10, min_samples=3)
    expression = median_ratio_log_expression(counts)

    samples = parse_soft_samples(RAW / "GSE190971_family.soft.gz")
    rows: list[dict[str, object]] = []
    for sample in samples:
        title = str(sample.get("title", ""))
        if not title.startswith("mRNA_") or not title.endswith("_PLAC"):
            continue
        specimen = title.removeprefix("mRNA_")
        disease = str(sample.get("disease_state_normal_or_pe", "")).upper()
        group = "PE" if disease == "PE" else "Control" if disease == "NP" else "Unclear"
        matrix_column = f"{specimen}_{'PE' if group == 'PE' else 'NORMAL'}"
        ga = parse_gestational_weeks(sample.get("gestational_age_at_delivery"))
        rows.append(
            {
                "cohort_id": "GSE190971",
                "original_study": "GSE190971",
                "geo_accession": sample["geo_accession"],
                "matrix_column": matrix_column,
                "group": group,
                "pe": 1 if group == "PE" else 0,
                "gestational_age_weeks": ga,
                "gestational_age_category": "Preterm" if np.isfinite(ga) and ga < 37 else "Term" if np.isfinite(ga) else "Unknown",
                "fetal_sex": clean_sex(sample.get("newborngender")),
                "chronic_hypertension": np.nan,
                "tissue": "placenta",
                "platform": "Illumina HiSeq 2000 RNA-seq (GPL11154)",
                "role": "external_replication",
                "independently_collected": True,
                "included_primary": matrix_column in expression.columns,
                "exclusion_reason": "" if matrix_column in expression.columns else "missing_from_placenta_count_matrix",
                "source_file": "GSE190971_Raw_gene_counts_matrix_PLAC.txt.gz",
                "provenance_notes": "Independent Oxford placenta RNA-seq cohort",
            }
        )
    manifest = pd.DataFrame(rows)
    columns = manifest.loc[manifest["included_primary"], "matrix_column"].tolist()
    return manifest, expression.loc[:, columns]


def read_ncbi_gene_map() -> dict[str, str]:
    columns = [
        "tax_id",
        "GeneID",
        "Symbol",
        "LocusTag",
        "Synonyms",
        "dbXrefs",
        "chromosome",
        "map_location",
        "description",
        "type_of_gene",
        "Symbol_from_nomenclature_authority",
        "Full_name_from_nomenclature_authority",
        "Nomenclature_status",
        "Other_designations",
        "Modification_date",
        "Feature_type",
    ]
    table = pd.read_csv(
        RAW / "Homo_sapiens.gene_info.gz",
        sep="\t",
        names=columns,
        header=0,
        dtype=str,
        low_memory=False,
    )
    return dict(zip(table["GeneID"], table["Symbol"].str.upper()))


def prepare_gse204835() -> tuple[pd.DataFrame, pd.DataFrame]:
    raw = pd.read_csv(RAW / "GSE204835_counts.csv.gz", dtype={"gene": str})
    mapping = read_ncbi_gene_map()
    raw["symbol"] = raw["gene"].map(mapping)
    mapping_rate = raw["symbol"].notna().mean()
    if mapping_rate < 0.70:
        raise RuntimeError(f"Unexpectedly low Entrez-to-symbol mapping rate: {mapping_rate:.3f}")
    counts = raw.dropna(subset=["symbol"]).drop(columns=["gene"]).groupby("symbol", sort=False).sum()
    counts = filter_count_matrix(counts, min_count=10, min_samples=3)
    expression = median_ratio_log_expression(counts)

    samples = parse_soft_samples(RAW / "GSE204835_family.soft.gz")
    rows: list[dict[str, object]] = []
    for sample in samples:
        matrix_column = first_description_token(sample.get("description"))
        disease = str(sample.get("disease_state", ""))
        ga_category = str(sample.get("gestational_age", "Unknown"))
        if disease == "Control":
            group = "Control"
        elif disease == "Preeclampsia":
            group = "PE"
        elif disease == "Chronic_hypertension":
            group = "Chronic_hypertension"
        elif disease == "Supermposed_preeclampsia":
            group = "Superimposed_PE"
        else:
            group = "Unclear"
        primary = matrix_column in expression.columns and (
            group == "Control" or (group == "PE" and ga_category == "Term")
        )
        reason = ""
        if matrix_column not in expression.columns:
            reason = "missing_from_count_matrix"
        elif not primary:
            reason = "not_in_term_de_novo_PE_primary_contrast"
        rows.append(
            {
                "cohort_id": "GSE204835",
                "original_study": "GSE204835",
                "geo_accession": sample["geo_accession"],
                "matrix_column": matrix_column,
                "group": group,
                "pe": 1 if group in {"PE", "Superimposed_PE"} else 0 if group in {"Control", "Chronic_hypertension"} else np.nan,
                "gestational_age_weeks": np.nan,
                "gestational_age_category": ga_category,
                "fetal_sex": "Unknown",
                "chronic_hypertension": int(group in {"Chronic_hypertension", "Superimposed_PE"}),
                "tissue": "FFPE placenta",
                "platform": "NextSeq 2000 QuantSeq 3-prime RNA-seq (GPL30173)",
                "role": "external_replication",
                "independently_collected": True,
                "included_primary": primary,
                "exclusion_reason": reason,
                "source_file": "GSE204835_counts.csv.gz",
                "provenance_notes": "Independent Michigan FFPE placenta cohort; primary contrast restricted to term de novo PE versus term controls",
            }
        )
    manifest = pd.DataFrame(rows)
    columns = [c for c in manifest["matrix_column"].dropna() if c in expression.columns]
    return manifest, expression.loc[:, columns]


def prepare_gse177049() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    genes = pd.read_excel(RAW / "GSE177049_expression.xlsx", sheet_name="gene expression", header=16)
    genes["Gene_Name"] = genes["Gene_Name"].astype(str).str.upper()
    gene_samples = [f"N{i}" for i in range(1, 6)] + [f"PE{i}" for i in range(1, 6)]
    gene_expr = genes.groupby("Gene_Name", sort=False)[gene_samples].mean()
    gene_expr = np.log2(gene_expr.astype(float) + 0.5)

    mirna = pd.read_excel(RAW / "GSE177049_expression.xlsx", sheet_name="miRNA expression", header=29)
    mirna["Mature_ID"] = mirna["Mature_ID"].astype(str).str.lower()
    mirna_expr = mirna.groupby("Mature_ID", sort=False)[gene_samples].mean()
    mirna_expr = np.log2(mirna_expr.astype(float) + 0.5)

    samples = parse_soft_samples(RAW / "GSE177049_family.soft.gz")
    rows: list[dict[str, object]] = []
    for sample in samples:
        title = str(sample.get("title", ""))
        if " - RNA Sequencing" not in title or "miRNA" in title:
            continue
        patient = title.split(" - ", 1)[0]
        group = "PE" if patient.startswith("PE") else "Control"
        rows.append(
            {
                "cohort_id": "GSE177049",
                "original_study": "GSE177049",
                "geo_accession": sample["geo_accession"],
                "matrix_column": patient,
                "group": group,
                "pe": 1 if group == "PE" else 0,
                "gestational_age_weeks": np.nan,
                "gestational_age_category": "Preterm",
                "fetal_sex": "Unknown",
                "chronic_hypertension": np.nan,
                "tissue": "placenta",
                "platform": "Illumina NovaSeq 6000 mRNA; paired NextSeq 500 miRNA",
                "role": "exploratory_paired_mrna_mirna",
                "independently_collected": True,
                "included_primary": False,
                "exclusion_reason": "exploratory_small_n_only",
                "source_file": "GSE177049_expression.xlsx",
                "provenance_notes": "Five early-onset PE and five normotensive preterm controls with paired mRNA/miRNA profiling",
            }
        )
    return pd.DataFrame(rows), gene_expr, mirna_expr


def save_matrix(matrix: pd.DataFrame, name: str) -> None:
    matrix.to_csv(PROCESSED / f"{name}.csv.gz", compression="gzip")


def main() -> None:
    new_manifest, historical_manifest, gse75010 = prepare_gse75010()
    external_190_manifest, external_190_expr = prepare_gse190971()
    external_204_manifest, external_204_expr = prepare_gse204835()
    exploratory_manifest, exploratory_gene, exploratory_mirna = prepare_gse177049()

    discovery_expr = gse75010["discovery"]
    discovery_expr.columns = discovery_expr.columns.get_level_values(0)
    historical_expr = gse75010["historical"]
    historical_expr.columns = historical_expr.columns.get_level_values(0)

    save_matrix(discovery_expr, "GSE75010_BioBank157_expression")
    save_matrix(historical_expr, "GSE75010_historical173_expression")
    save_matrix(external_190_expr, "GSE190971_expression")
    save_matrix(external_204_expr, "GSE204835_expression")
    save_matrix(exploratory_gene, "GSE177049_gene_expression")
    save_matrix(exploratory_mirna, "GSE177049_mirna_expression")

    manifest = pd.concat(
        [new_manifest, historical_manifest, external_190_manifest, external_204_manifest, exploratory_manifest],
        ignore_index=True,
    )
    manifest.to_csv(PROCESSED / "cohort_manifest.csv", index=False)

    summary = (
        manifest.groupby(["cohort_id", "role", "group"], dropna=False)
        .agg(n_records=("matrix_column", "size"), n_primary=("included_primary", "sum"))
        .reset_index()
    )
    summary.to_csv(PROCESSED / "cohort_counts.csv", index=False)

    raw_files = sorted(path for path in RAW.iterdir() if path.is_file() and not path.name.startswith("."))
    file_manifest = {
        "files": [
            {"path": str(path.relative_to(ROOT)).replace("\\", "/"), "bytes": path.stat().st_size, "sha256": sha256_file(path)}
            for path in raw_files
        ],
        "cohort_counts": summary.to_dict(orient="records"),
        "notes": {
            "GSE75010": "The integrated matrix has 330 columns; only the 157 RCWIH BioBank samples are discovery data.",
            "historical": "The seven imported studies sum to 173 selected PE/control samples and are not external validation.",
            "GSE204835_primary": "Term de novo PE versus term normotensive controls; other hypertensive groups are sensitivity-only.",
            "GSE204835_metadata_discrepancy": "The GEO series summary states 11 term and 5 preterm de novo PE, while current sample-level SOFT fields label 12 term and 4 preterm; this workflow follows sample-level fields and records the discrepancy.",
        },
    }
    write_json(file_manifest, PROCESSED / "data_manifest.json")

    print(summary.to_string(index=False))
    print("\nProcessed matrices:")
    for path in sorted(PROCESSED.glob("*_expression.csv.gz")):
        print(f"  {path.name}: {path.stat().st_size:,} bytes")


if __name__ == "__main__":
    main()
