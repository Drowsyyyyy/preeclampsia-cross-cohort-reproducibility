"""Run gene-, pathway-, meta-, sensitivity-, and negative-control analyses."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import stats

from common import (
    bh_fdr,
    dersimonian_laird,
    fit_matrix_ols,
    read_gmt,
    score_gene_sets,
    write_json,
    zscore_rows,
)


ROOT = Path(__file__).resolve().parents[1]
PROCESSED = ROOT / "data" / "processed"
RAW = ROOT / "data" / "raw"
RESULTS = ROOT / "results"
TABLES = RESULTS / "tables"
COHORT_RESULTS = RESULTS / "cohort_level_results"
TABLES.mkdir(parents=True, exist_ok=True)
COHORT_RESULTS.mkdir(parents=True, exist_ok=True)


def load_plan() -> dict:
    return yaml.safe_load((ROOT / "config" / "analysis_plan.yaml").read_text(encoding="utf-8"))


def load_expression(name: str) -> pd.DataFrame:
    expression = pd.read_csv(PROCESSED / f"{name}_expression.csv.gz", index_col=0)
    expression.index = expression.index.astype(str).str.upper()
    return expression.groupby(expression.index, sort=False).mean()


def metadata_for(manifest: pd.DataFrame, cohort: str, primary_only: bool = False) -> pd.DataFrame:
    meta = manifest.loc[manifest["cohort_id"].eq(cohort)].copy()
    if primary_only:
        meta = meta.loc[meta["included_primary"].astype(bool)]
    return meta.set_index("matrix_column", drop=False)


def discovery_design(meta: pd.DataFrame) -> pd.DataFrame:
    design = pd.DataFrame(index=meta.index)
    design["intercept"] = 1.0
    design["pe"] = meta["pe"].astype(float)
    design["gestational_age_centered"] = meta["gestational_age_weeks"].astype(float) - meta["gestational_age_weeks"].astype(float).mean()
    design["male"] = meta["fetal_sex"].eq("M").astype(float)
    design["chronic_hypertension"] = meta["chronic_hypertension"].astype(float)
    return design


def group_design(meta: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame({"intercept": 1.0, "pe": meta["pe"].astype(float)}, index=meta.index)


def group_ga_design(meta: pd.DataFrame) -> pd.DataFrame:
    ga = meta["gestational_age_weeks"].astype(float)
    keep = ga.notna() & meta["pe"].notna()
    design = pd.DataFrame(index=meta.index[keep])
    design["intercept"] = 1.0
    design["pe"] = meta.loc[keep, "pe"].astype(float)
    design["gestational_age_centered"] = ga.loc[keep] - ga.loc[keep].mean()
    return design


def selected_gene_sets(plan: dict) -> dict[str, list[str]]:
    hallmark = read_gmt(RAW / "MSigDB_Hallmark_2020.gmt")
    chosen: dict[str, list[str]] = {}
    for name in plan["primary_pathways"]:
        if name not in hallmark:
            raise KeyError(f"Missing prespecified hallmark gene set: {name}")
        chosen[f"HALLMARK_{name.upper().replace(' ', '_').replace('-', '_')}"] = hallmark[name]
    for name, genes in plan["curated_gene_sets"].items():
        chosen[f"CURATED_{name}"] = genes
    return chosen


def run_cohort(
    cohort_name: str,
    expression: pd.DataFrame,
    design: pd.DataFrame,
    gene_sets: dict[str, list[str]],
) -> dict[str, pd.DataFrame]:
    samples = [sample for sample in expression.columns if sample in design.index]
    expression = expression.loc[:, samples]
    design = design.loc[samples]

    standardized = zscore_rows(expression).dropna(how="all")
    gene_results = fit_matrix_ols(standardized, design, "pe")
    gene_results.insert(0, "cohort", cohort_name)

    pathway_scores, pathway_audit = score_gene_sets(expression, gene_sets)
    pathway_scores = zscore_rows(pathway_scores).dropna(how="all")
    pathway_results = fit_matrix_ols(pathway_scores, design, "pe")
    pathway_results.insert(0, "cohort", cohort_name)
    pathway_audit.insert(0, "cohort", cohort_name)

    sample_scores = pathway_scores.T.copy()
    sample_scores.index.name = "matrix_column"
    sample_scores.reset_index().to_csv(COHORT_RESULTS / f"{cohort_name}_pathway_scores.csv", index=False)
    gene_results.to_csv(COHORT_RESULTS / f"{cohort_name}_gene_effects.csv", index=False)
    pathway_results.to_csv(COHORT_RESULTS / f"{cohort_name}_pathway_effects.csv", index=False)
    pathway_audit.to_csv(COHORT_RESULTS / f"{cohort_name}_gene_set_coverage.csv", index=False)
    return {
        "genes": gene_results,
        "pathways": pathway_results,
        "coverage": pathway_audit,
        "scores": pathway_scores,
    }


def meta_analyze(results: list[pd.DataFrame], feature_type: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    stacked = pd.concat(results, ignore_index=True)
    rows: list[dict[str, object]] = []
    loco_rows: list[dict[str, object]] = []
    for feature, group in stacked.groupby("feature", sort=False):
        meta = dersimonian_laird(group["estimate"].to_numpy(), group["std_error"].to_numpy())
        rows.append({"feature_type": feature_type, "feature": feature, **meta})
        for omitted in group["cohort"].unique():
            retained = group.loc[~group["cohort"].eq(omitted)]
            loco = dersimonian_laird(retained["estimate"].to_numpy(), retained["std_error"].to_numpy())
            loco_rows.append(
                {
                    "feature_type": feature_type,
                    "feature": feature,
                    "omitted_cohort": omitted,
                    **loco,
                }
            )
    meta_df = pd.DataFrame(rows)
    meta_df["fdr"] = bh_fdr(meta_df["p_value"])
    meta_df = meta_df.sort_values("p_value", kind="stable").reset_index(drop=True)
    return meta_df, pd.DataFrame(loco_rows)


def matched_random_set(
    match_table: pd.DataFrame,
    pools: dict[tuple[int, int], np.ndarray],
    genes: list[str],
    rng: np.random.Generator,
) -> list[str]:
    requested = [g.upper() for g in genes if g.upper() in match_table.index]
    selected: list[str] = []
    excluded = set(requested)
    selected_set: set[str] = set()
    all_genes = match_table.index.to_numpy()
    for gene in requested:
        row = match_table.loc[gene]
        key = (int(row["mean_bin"]), int(row["sd_bin"]))
        pool = pools.get(key, all_genes)
        choice: str | None = None
        for _ in range(50):
            candidate = str(rng.choice(pool))
            if candidate not in excluded and candidate not in selected_set:
                choice = candidate
                break
        if choice is None:
            valid = [g for g in pool if g not in excluded and g not in selected_set]
            if not valid:
                valid = [g for g in all_genes if g not in excluded and g not in selected_set]
            choice = str(rng.choice(valid))
        selected.append(choice)
        selected_set.add(choice)
    return selected


def negative_controls(
    expression: pd.DataFrame,
    design: pd.DataFrame,
    gene_sets: dict[str, list[str]],
    observed: pd.DataFrame,
    n_random: int,
    n_permutations: int,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(seed)
    observed_map = observed.set_index("feature")["estimate"].to_dict()
    random_rows: list[dict[str, object]] = []
    permutation_rows: list[dict[str, object]] = []
    samples = [sample for sample in expression.columns if sample in design.index]
    expression = expression.loc[:, samples]
    design = design.loc[samples]
    expression.index = expression.index.astype(str).str.upper()
    expression = expression.groupby(expression.index, sort=False).mean()
    z_expression = zscore_rows(expression).dropna(how="all")
    match_table = pd.DataFrame(
        {"mean": expression.mean(axis=1), "sd": expression.std(axis=1, ddof=1)},
        index=expression.index,
    ).replace([np.inf, -np.inf], np.nan).dropna()
    match_table["mean_bin"] = pd.qcut(match_table["mean"], q=10, labels=False, duplicates="drop")
    match_table["sd_bin"] = pd.qcut(match_table["sd"], q=10, labels=False, duplicates="drop")
    pools = {
        (int(mean_bin), int(sd_bin)): group.index.to_numpy()
        for (mean_bin, sd_bin), group in match_table.groupby(["mean_bin", "sd_bin"], observed=True)
    }

    for name, genes in gene_sets.items():
        for iteration in range(n_random):
            sampled = matched_random_set(match_table, pools, genes, rng)
            if len(sampled) < 3:
                continue
            score = z_expression.loc[sampled].mean(axis=0).to_frame(name="random").T
            score = zscore_rows(score).dropna(how="all")
            effect = fit_matrix_ols(score, design, "pe").iloc[0]
            random_rows.append(
                {
                    "gene_set": name,
                    "iteration": iteration,
                    "estimate": effect["estimate"],
                    "n_genes": len(sampled),
                }
            )

        present = [gene.upper() for gene in genes if gene.upper() in z_expression.index]
        score = z_expression.loc[present].mean(axis=0).to_frame(name=name).T
        score = zscore_rows(score).dropna(how="all")
        for iteration in range(n_permutations):
            permuted = design.copy()
            permuted["pe"] = rng.permutation(permuted["pe"].to_numpy())
            effect = fit_matrix_ols(score, permuted, "pe").iloc[0]
            permutation_rows.append(
                {"gene_set": name, "iteration": iteration, "estimate": effect["estimate"]}
            )

    random_df = pd.DataFrame(random_rows)
    permutation_df = pd.DataFrame(permutation_rows)
    summaries: list[dict[str, object]] = []
    for name, observed_estimate in observed_map.items():
        random_null = random_df.loc[random_df["gene_set"].eq(name), "estimate"].dropna().to_numpy()
        perm_null = permutation_df.loc[permutation_df["gene_set"].eq(name), "estimate"].dropna().to_numpy()
        summaries.append(
            {
                "gene_set": name,
                "observed_estimate": observed_estimate,
                "matched_random_plus_one_p": (1 + np.sum(np.abs(random_null) >= abs(observed_estimate))) / (len(random_null) + 1),
                "label_permutation_plus_one_p": (1 + np.sum(np.abs(perm_null) >= abs(observed_estimate))) / (len(perm_null) + 1),
                "n_matched_random": len(random_null),
                "n_label_permutations": len(perm_null),
            }
        )
    return pd.DataFrame(summaries), pd.concat(
        {
            "matched_random": random_df.set_index(["gene_set", "iteration"]),
            "label_permutation": permutation_df.set_index(["gene_set", "iteration"]),
        }
    ).reset_index()


def composition_sensitivity(
    expression: pd.DataFrame,
    meta: pd.DataFrame,
    plan: dict,
    primary_gene_sets: dict[str, list[str]],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    proxy_scores, proxy_audit = score_gene_sets(expression, plan["cell_composition_proxies"])
    design = discovery_design(meta)
    for proxy in ["ENDOTHELIAL_MARKERS", "GENERAL_TROPHOBLAST_MARKERS"]:
        design[f"proxy_{proxy.lower()}"] = proxy_scores.loc[proxy, design.index]

    target_genes = [g for g in plan["target_genes"] if g in expression.index]
    target_expression = zscore_rows(expression.loc[target_genes]).dropna(how="all")
    target_results = fit_matrix_ols(target_expression, design, "pe")
    target_results.insert(0, "model", "composition_proxy_adjusted")

    pathway_scores, _ = score_gene_sets(expression, primary_gene_sets)
    pathway_scores = zscore_rows(pathway_scores).dropna(how="all")
    pathway_results = fit_matrix_ols(pathway_scores, design, "pe")
    pathway_results.insert(0, "model", "composition_proxy_adjusted")

    sample_proxy = proxy_scores.T.copy()
    sample_proxy.index.name = "matrix_column"
    sample_proxy = sample_proxy.reset_index()
    return target_results, pathway_results, sample_proxy.merge(meta.reset_index(drop=True), on="matrix_column", how="left"), proxy_audit


def exploratory_mirna_analysis(plan: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    manifest = pd.read_csv(PROCESSED / "cohort_manifest.csv")
    meta = metadata_for(manifest, "GSE177049", primary_only=False)
    gene = load_expression("GSE177049_gene")
    mirna = pd.read_csv(PROCESSED / "GSE177049_mirna_expression.csv.gz", index_col=0)
    mirna.index = mirna.index.astype(str).str.lower()
    design = group_design(meta)

    selected_mirnas = ["hsa-mir-31-5p", "hsa-mir-155-5p", "hsa-mir-214-3p", "hsa-mir-1290"]
    available_mirnas = [m for m in selected_mirnas if m in mirna.index]
    mirna_results = fit_matrix_ols(zscore_rows(mirna.loc[available_mirnas]).dropna(how="all"), design, "pe")
    mirna_results.insert(0, "analysis", "GSE177049_standardized_group_effect")

    pairs = [
        ("hsa-mir-214-3p", "PGF"),
        ("hsa-mir-214-3p", "NOS3"),
        ("hsa-mir-31-5p", "NOS3"),
        ("hsa-mir-155-5p", "NOS3"),
    ]
    correlation_rows: list[dict[str, object]] = []
    for mir, target in pairs:
        if mir not in mirna.index or target not in gene.index:
            correlation_rows.append({"mirna": mir, "target": target, "rho": np.nan, "p_value": np.nan, "n": 0})
            continue
        samples = [sample for sample in gene.columns if sample in mirna.columns]
        rho, p = stats.spearmanr(mirna.loc[mir, samples], gene.loc[target, samples])
        correlation_rows.append({"mirna": mir, "target": target, "rho": rho, "p_value": p, "n": len(samples)})
    correlations = pd.DataFrame(correlation_rows)
    correlations["fdr"] = bh_fdr(correlations["p_value"])
    return mirna_results, correlations


def main() -> None:
    plan = load_plan()
    manifest = pd.read_csv(PROCESSED / "cohort_manifest.csv")
    gene_sets = selected_gene_sets(plan)

    discovery_expr = load_expression("GSE75010_BioBank157")
    discovery_meta = metadata_for(manifest, "GSE75010_BioBank157", primary_only=True)
    disc_design = discovery_design(discovery_meta)

    # Full adjusted discovery differential-expression table in original log2 units.
    discovery_de = fit_matrix_ols(discovery_expr, disc_design, "pe")
    discovery_de.insert(0, "model", "adjusted_primary")
    discovery_de.to_csv(TABLES / "discovery_gene_differential_expression.csv.gz", index=False, compression="gzip")

    independent: dict[str, dict[str, pd.DataFrame]] = {}
    independent["GSE75010_BioBank157"] = run_cohort(
        "GSE75010_BioBank157", discovery_expr, disc_design, gene_sets
    )

    expr_190 = load_expression("GSE190971")
    meta_190 = metadata_for(manifest, "GSE190971", primary_only=True)
    independent["GSE190971"] = run_cohort("GSE190971", expr_190, group_design(meta_190), gene_sets)

    expr_204 = load_expression("GSE204835")
    meta_204 = metadata_for(manifest, "GSE204835", primary_only=True)
    independent["GSE204835_term"] = run_cohort(
        "GSE204835_term", expr_204, group_design(meta_204), gene_sets
    )

    # External sensitivity models.
    if meta_190["gestational_age_weeks"].notna().sum() >= 10:
        ga_design_190 = group_ga_design(meta_190)
        run_cohort("GSE190971_GA_adjusted", expr_190, ga_design_190, gene_sets)

    meta_204_all = metadata_for(manifest, "GSE204835", primary_only=False)
    meta_204_all = meta_204_all.loc[meta_204_all["group"].isin(["PE", "Control"])]
    design_204_all = group_design(meta_204_all)
    design_204_all["preterm"] = meta_204_all["gestational_age_category"].eq("Preterm").astype(float)
    run_cohort("GSE204835_all_de_novo", expr_204, design_204_all, gene_sets)

    # Historical cohorts are analyzed separately and never enter the external meta-analysis.
    historical_expr = load_expression("GSE75010_historical173")
    historical_gene_results: list[pd.DataFrame] = []
    historical_path_results: list[pd.DataFrame] = []
    for cohort in plan["cohorts"]["historical_internal_replication"]["ids"]:
        meta = metadata_for(manifest, cohort, primary_only=True)
        expr = historical_expr.loc[:, [c for c in meta.index if c in historical_expr.columns]]
        ga_available = meta["gestational_age_weeks"].notna().mean() >= 0.80
        design = group_ga_design(meta) if ga_available else group_design(meta)
        result = run_cohort(cohort, expr, design, gene_sets)
        historical_gene_results.append(result["genes"])
        historical_path_results.append(result["pathways"])

    historical_genes = pd.concat(historical_gene_results, ignore_index=True)
    historical_paths = pd.concat(historical_path_results, ignore_index=True)
    historical_genes.to_csv(TABLES / "historical_internal_gene_effects.csv.gz", index=False, compression="gzip")
    historical_paths.to_csv(TABLES / "historical_internal_pathway_effects.csv", index=False)

    # Random-effects meta-analysis includes only independently collected resources.
    gene_meta, gene_loco = meta_analyze([x["genes"] for x in independent.values()], "gene")
    pathway_meta, pathway_loco = meta_analyze([x["pathways"] for x in independent.values()], "pathway")
    gene_meta.to_csv(TABLES / "independent_gene_meta_analysis.csv.gz", index=False, compression="gzip")
    pathway_meta.to_csv(TABLES / "independent_pathway_meta_analysis.csv", index=False)
    pd.concat([gene_loco, pathway_loco], ignore_index=True).to_csv(TABLES / "leave_one_cohort_out_meta_analysis.csv.gz", index=False, compression="gzip")

    # Prespecified target-gene summary.
    targets = {gene.upper() for gene in plan["target_genes"]}
    target_cohort = pd.concat([x["genes"] for x in independent.values()], ignore_index=True)
    target_cohort = target_cohort.loc[target_cohort["feature"].isin(targets)]
    target_meta = gene_meta.loc[gene_meta["feature"].isin(targets)]
    target_cohort.to_csv(TABLES / "target_gene_cohort_effects.csv", index=False)
    target_meta.to_csv(TABLES / "target_gene_meta_analysis.csv", index=False)

    # Composition-proxy sensitivity analysis in discovery data.
    target_sens, path_sens, proxy_samples, proxy_audit = composition_sensitivity(
        discovery_expr, discovery_meta, plan, gene_sets
    )
    target_sens.to_csv(TABLES / "composition_adjusted_target_genes.csv", index=False)
    path_sens.to_csv(TABLES / "composition_adjusted_pathways.csv", index=False)
    proxy_samples.to_csv(TABLES / "discovery_cell_composition_proxy_scores.csv", index=False)
    proxy_audit.to_csv(TABLES / "cell_composition_proxy_coverage.csv", index=False)

    # Discovery negative controls.
    neg_cfg = plan["statistics"]["negative_controls"]
    negative_summary, negative_full = negative_controls(
        discovery_expr,
        disc_design,
        gene_sets,
        independent["GSE75010_BioBank157"]["pathways"],
        int(neg_cfg["matched_random_gene_sets"]),
        int(neg_cfg["outcome_label_permutations"]),
        int(neg_cfg["random_seed"]),
    )
    negative_summary.to_csv(TABLES / "negative_control_summary.csv", index=False)
    negative_full.to_csv(TABLES / "negative_control_iterations.csv.gz", index=False, compression="gzip")

    mirna_results, mirna_correlations = exploratory_mirna_analysis(plan)
    mirna_results.to_csv(TABLES / "GSE177049_mirna_effects.csv", index=False)
    mirna_correlations.to_csv(TABLES / "GSE177049_mirna_target_correlations.csv", index=False)

    # Compact, machine-readable run summary.
    summary = {
        "independent_cohorts": {
            name: {
                "n_samples": int(result["genes"]["n"].iloc[0]),
                "n_genes": int(len(result["genes"])),
                "n_gene_sets": int(len(result["pathways"])),
            }
            for name, result in independent.items()
        },
        "discovery_significant_genes_fdr_0_05": int((discovery_de["fdr"] < 0.05).sum()),
        "independent_meta_significant_genes_fdr_0_05": int((gene_meta["fdr"] < 0.05).sum()),
        "independent_meta_significant_pathways_fdr_0_05": int((pathway_meta["fdr"] < 0.05).sum()),
        "negative_control_iterations_per_set": {
            "matched_random": int(neg_cfg["matched_random_gene_sets"]),
            "label_permutation": int(neg_cfg["outcome_label_permutations"]),
        },
    }
    write_json(summary, RESULTS / "analysis_summary.json")

    print(json.dumps(summary, indent=2))
    print("\nIndependent pathway meta-analysis:")
    print(pathway_meta.to_string(index=False))
    print("\nTarget-gene meta-analysis:")
    print(target_meta.sort_values("p_value").to_string(index=False))


if __name__ == "__main__":
    main()
