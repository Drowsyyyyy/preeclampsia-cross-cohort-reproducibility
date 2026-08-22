"""Shared utilities for the preeclampsia cross-cohort study.

The project intentionally avoids hidden state: parsing, normalization, model fitting,
multiple-testing correction, and meta-analysis are implemented here and called by
the numbered scripts.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import math
import re
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy import stats


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(obj: object, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, sort_keys=True), encoding="utf-8")


def parse_soft_samples(path: Path) -> list[dict[str, object]]:
    """Parse sample-level fields from a compressed GEO family SOFT file."""
    samples: list[dict[str, object]] = []
    current: dict[str, object] | None = None
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\r\n")
            if line.startswith("^SAMPLE = "):
                current = {"_sample_accession": line.split(" = ", 1)[1]}
                samples.append(current)
                continue
            if current is None or not line.startswith("!Sample_"):
                continue
            key, value = line.split(" = ", 1)
            key = key.removeprefix("!Sample_")
            current.setdefault(key, [])
            assert isinstance(current[key], list)
            current[key].append(value)

    flattened: list[dict[str, object]] = []
    for sample in samples:
        out: dict[str, object] = {"geo_accession": sample["_sample_accession"]}
        for key, values in sample.items():
            if key == "_sample_accession":
                continue
            assert isinstance(values, list)
            if key == "characteristics_ch1":
                for item in values:
                    if ": " in item:
                        char_key, char_value = item.split(": ", 1)
                        out[normalize_key(char_key)] = char_value.strip()
                out["characteristics_raw"] = " | ".join(values)
            elif len(values) == 1:
                out[key] = values[0]
            else:
                out[key] = " | ".join(values)
        flattened.append(out)
    return flattened


def normalize_key(value: str) -> str:
    value = value.strip().lower()
    value = re.sub(r"[^a-z0-9]+", "_", value)
    return value.strip("_")


def first_description_token(value: object) -> str | None:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return None
    return str(value).split(" | ", 1)[0].strip()


def parse_number(value: object) -> float:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return np.nan
    text = str(value).strip().lower()
    if text in {"", "na", "n/a", "none", "nan", "unknown", "not available"}:
        return np.nan
    match = re.search(r"[-+]?\d+(?:\.\d+)?", text)
    return float(match.group()) if match else np.nan


def parse_gestational_weeks(value: object) -> float:
    """Convert strings such as '38', '41+2', or '35 weeks 3 days' to weeks."""
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return np.nan
    text = str(value).strip().lower()
    if text in {"", "na", "n/a", "none", "nan", "term", "preterm"}:
        return np.nan
    plus = re.fullmatch(r"\s*(\d+)\s*\+\s*(\d+)\s*", text)
    if plus:
        return float(plus.group(1)) + float(plus.group(2)) / 7.0
    weeks_days = re.search(r"(\d+(?:\.\d+)?)\s*(?:weeks?|w)", text)
    if weeks_days:
        weeks = float(weeks_days.group(1))
        days = re.search(r"(\d+(?:\.\d+)?)\s*(?:days?|d)", text)
        return weeks + (float(days.group(1)) / 7.0 if days else 0.0)
    return parse_number(text)


def bh_fdr(p_values: Iterable[float]) -> np.ndarray:
    p = np.asarray(list(p_values), dtype=float)
    out = np.full(p.shape, np.nan)
    valid = np.isfinite(p)
    if not valid.any():
        return out
    pv = p[valid]
    order = np.argsort(pv)
    ranked = pv[order]
    n = len(ranked)
    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    restored = np.empty_like(adjusted)
    restored[order] = adjusted
    out[valid] = restored
    return out


def median_ratio_log_expression(counts: pd.DataFrame) -> pd.DataFrame:
    """DESeq-style median-ratio size factors followed by log2(count/sf + 0.5)."""
    numeric = counts.astype(float)
    # The canonical median-ratio estimator uses genes positive in every sample.
    # Treating structural zeros as missing would bias sparse matrices toward a
    # zero size factor, so we use all-positive genes and fall back to library size.
    usable = (numeric > 0).all(axis=1)
    if usable.sum() >= 100:
        geometric_means = np.exp(np.log(numeric.loc[usable]).mean(axis=1))
        ratios = numeric.loc[usable].div(geometric_means, axis=0)
        size_factors = ratios.median(axis=0, skipna=True)
    else:
        library_sizes = numeric.sum(axis=0)
        size_factors = library_sizes / np.exp(np.mean(np.log(library_sizes)))
    if (~np.isfinite(size_factors)).any() or (size_factors <= 0).any():
        library_sizes = numeric.sum(axis=0)
        size_factors = library_sizes / np.exp(np.mean(np.log(library_sizes)))
    size_factors = size_factors / np.exp(np.mean(np.log(size_factors)))
    normalized = numeric.div(size_factors, axis=1)
    return np.log2(normalized + 0.5)


def filter_count_matrix(counts: pd.DataFrame, min_count: int = 10, min_samples: int = 3) -> pd.DataFrame:
    keep = (counts >= min_count).sum(axis=1) >= min_samples
    return counts.loc[keep].copy()


def fit_matrix_ols(
    expression: pd.DataFrame,
    design: pd.DataFrame,
    coefficient: str,
) -> pd.DataFrame:
    """Fit ordinary least squares to all expression rows using matrix algebra.

    Expression is genes/features x samples; design is samples x covariates.
    Classical standard errors are reported. The function requires complete data.
    """
    common = [sample for sample in expression.columns if sample in design.index]
    if not common:
        raise ValueError("No shared samples between expression and design")
    xdf = design.loc[common].astype(float)
    if xdf.isna().any().any():
        raise ValueError("Design matrix contains missing values")
    ydf = expression.loc[:, common].astype(float)
    if ydf.isna().any().any():
        ydf = ydf.apply(lambda row: row.fillna(row.median()), axis=1)
    x = xdf.to_numpy()
    y = ydf.to_numpy().T
    xtx_inv = np.linalg.pinv(x.T @ x)
    beta = xtx_inv @ x.T @ y
    resid = y - x @ beta
    rank = np.linalg.matrix_rank(x)
    df_resid = x.shape[0] - rank
    if df_resid <= 0:
        raise ValueError("Non-positive residual degrees of freedom")
    sigma2 = np.sum(resid**2, axis=0) / df_resid
    coef_index = list(xdf.columns).index(coefficient)
    se = np.sqrt(np.maximum(sigma2 * xtx_inv[coef_index, coef_index], 0))
    estimate = beta[coef_index]
    with np.errstate(divide="ignore", invalid="ignore"):
        t_stat = estimate / se
    p_value = 2 * stats.t.sf(np.abs(t_stat), df=df_resid)
    result = pd.DataFrame(
        {
            "feature": ydf.index.astype(str),
            "estimate": estimate,
            "std_error": se,
            "t_stat": t_stat,
            "p_value": p_value,
            "n": x.shape[0],
            "df_resid": df_resid,
        }
    )
    result["fdr"] = bh_fdr(result["p_value"])
    return result.sort_values("p_value", kind="stable").reset_index(drop=True)


def zscore_rows(expression: pd.DataFrame) -> pd.DataFrame:
    means = expression.mean(axis=1)
    stds = expression.std(axis=1, ddof=1).replace(0, np.nan)
    return expression.sub(means, axis=0).div(stds, axis=0)


def score_gene_sets(expression: pd.DataFrame, gene_sets: dict[str, list[str]]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Mean within-cohort gene z-scores for each gene set and membership audit."""
    expr = expression.copy()
    expr.index = expr.index.astype(str).str.upper()
    expr = expr.groupby(expr.index, sort=False).mean()
    z = zscore_rows(expr)
    scores: dict[str, pd.Series] = {}
    audit: list[dict[str, object]] = []
    for name, genes in gene_sets.items():
        requested = sorted({str(g).upper() for g in genes})
        present = [g for g in requested if g in z.index and z.loc[g].notna().any()]
        if len(present) >= 3:
            scores[name] = z.loc[present].mean(axis=0)
        audit.append(
            {
                "gene_set": name,
                "n_requested": len(requested),
                "n_present": len(present),
                "present_genes": ";".join(present),
                "missing_genes": ";".join(sorted(set(requested) - set(present))),
            }
        )
    return pd.DataFrame(scores).T, pd.DataFrame(audit)


def dersimonian_laird(effects: np.ndarray, standard_errors: np.ndarray) -> dict[str, float]:
    effects = np.asarray(effects, dtype=float)
    standard_errors = np.asarray(standard_errors, dtype=float)
    valid = np.isfinite(effects) & np.isfinite(standard_errors) & (standard_errors > 0)
    y = effects[valid]
    se = standard_errors[valid]
    k = len(y)
    if k == 0:
        return {key: np.nan for key in ["estimate", "std_error", "z", "p_value", "ci_low", "ci_high", "tau2", "i2", "q", "k"]}
    w = 1.0 / se**2
    fixed = np.sum(w * y) / np.sum(w)
    q = np.sum(w * (y - fixed) ** 2)
    c = np.sum(w) - np.sum(w**2) / np.sum(w)
    tau2 = max(0.0, (q - (k - 1)) / c) if k > 1 and c > 0 else 0.0
    wr = 1.0 / (se**2 + tau2)
    estimate = np.sum(wr * y) / np.sum(wr)
    meta_se = math.sqrt(1.0 / np.sum(wr))
    z = estimate / meta_se
    p = 2 * stats.norm.sf(abs(z))
    i2 = max(0.0, (q - (k - 1)) / q) * 100 if k > 1 and q > 0 else 0.0
    return {
        "estimate": estimate,
        "std_error": meta_se,
        "z": z,
        "p_value": p,
        "ci_low": estimate - 1.96 * meta_se,
        "ci_high": estimate + 1.96 * meta_se,
        "tau2": tau2,
        "i2": i2,
        "q": q,
        "k": float(k),
    }


def read_gmt(path: Path) -> dict[str, list[str]]:
    gene_sets: dict[str, list[str]] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        parts = line.rstrip().split("\t")
        if len(parts) >= 3:
            gene_sets[parts[0]] = parts[2:]
    return gene_sets
