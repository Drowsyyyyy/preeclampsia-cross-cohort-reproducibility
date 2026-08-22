"""Create publication-style figures from retained analysis artifacts."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "results" / "tables"
COHORT = ROOT / "results" / "cohort_level_results"
FIGURES = ROOT / "results" / "figures"
FIGURES.mkdir(parents=True, exist_ok=True)

sns.set_theme(style="whitegrid", context="talk")
plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "axes.titleweight": "bold",
        "figure.dpi": 160,
        "savefig.dpi": 220,
        "savefig.bbox": "tight",
    }
)


def save(fig: plt.Figure, name: str) -> None:
    fig.savefig(FIGURES / name, facecolor="white")
    plt.close(fig)


def study_design() -> None:
    fig, ax = plt.subplots(figsize=(13, 7))
    ax.set_xlim(0, 13)
    ax.set_ylim(0, 7)
    ax.axis("off")
    ax.set_title("Study design and evidence hierarchy", fontsize=22, pad=18)

    boxes = [
        (0.4, 4.7, 3.5, 1.35, "Discovery\nGSE75010 new BioBank\n80 PE / 77 non-PE", "#d9edf7"),
        (4.7, 4.7, 3.5, 1.35, "External replication 1\nGSE190971 fresh placenta\n7 PE / 6 controls", "#dff0d8"),
        (9.0, 4.7, 3.5, 1.35, "External replication 2\nGSE204835 FFPE, term-only\n12 PE / 12 controls", "#dff0d8"),
        (0.4, 1.8, 5.0, 1.35, "Historical internal replication\n7 earlier GEO studies, 173 samples\nAlready embedded in GSE75010", "#fce5cd"),
        (7.6, 1.8, 4.9, 1.35, "Exploratory paired mRNA/miRNA\nGSE177049\n5 early-onset PE / 5 preterm controls", "#eadcf8"),
    ]
    for x, y, w, h, text, color in boxes:
        ax.add_patch(plt.Rectangle((x, y), w, h, facecolor=color, edgecolor="#455a64", linewidth=1.8))
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=12.5)

    ax.annotate("", xy=(4.55, 5.38), xytext=(3.95, 5.38), arrowprops=dict(arrowstyle="->", lw=2, color="#455a64"))
    ax.annotate("", xy=(8.85, 5.38), xytext=(8.25, 5.38), arrowprops=dict(arrowstyle="->", lw=2, color="#455a64"))
    ax.annotate("Independent random-effects meta-analysis", xy=(6.5, 4.25), ha="center", fontsize=12, color="#263238")
    ax.plot([2.9, 2.9], [3.15, 4.6], color="#9e9e9e", linestyle="--")
    ax.text(2.9, 3.55, "Not external validation", ha="center", va="center", fontsize=11, color="#8a4b08", bbox=dict(facecolor="white", edgecolor="none", pad=2))
    ax.text(10.05, 0.85, "miRNA findings treated as descriptive only", ha="center", fontsize=11, color="#5e3c8a")
    ax.text(6.5, 0.25, "All cohorts were processed separately; no raw-scale pooling across platforms", ha="center", fontsize=12, weight="bold")
    save(fig, "figure1_study_design.png")


def pathway_forest() -> None:
    cohorts = ["GSE75010_BioBank157", "GSE190971", "GSE204835_term"]
    labels = {"GSE75010_BioBank157": "Discovery", "GSE190971": "External 1", "GSE204835_term": "External 2"}
    frames = []
    for cohort in cohorts:
        frame = pd.read_csv(COHORT / f"{cohort}_pathway_effects.csv")
        frame["source"] = labels[cohort]
        frames.append(frame)
    data = pd.concat(frames, ignore_index=True)
    meta = pd.read_csv(TABLES / "independent_pathway_meta_analysis.csv")
    meta = meta.rename(columns={"std_error": "std_error"})
    meta["source"] = "Random-effects meta"
    data = pd.concat([data, meta[data.columns.intersection(meta.columns).tolist()]], ignore_index=True)

    order = pd.read_csv(TABLES / "independent_pathway_meta_analysis.csv").sort_values("estimate")["feature"].tolist()
    pretty = {
        "CURATED_ANTIANGIOGENIC_VASCULAR_STRESS": "Antiangiogenic / vascular stress",
        "HALLMARK_HYPOXIA": "Hypoxia",
        "HALLMARK_TNF_ALPHA_SIGNALING_VIA_NF_KB": "TNF-alpha / NF-kB",
        "HALLMARK_INFLAMMATORY_RESPONSE": "Inflammatory response",
        "CURATED_HIF1_TWIST1_MIR214_HOST_AXIS": "HIF1 / TWIST1 / miR-214 host",
        "CURATED_EXTRAVILLOUS_TROPHOBLAST_INVASION": "Extravillous trophoblast invasion",
        "CURATED_PROANGIOGENIC_SUPPORT": "Proangiogenic support",
        "CURATED_ENDOTHELIAL_NO_SIGNALING": "Endothelial NO signaling",
        "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION": "Epithelial-mesenchymal transition",
        "HALLMARK_ANGIOGENESIS": "Angiogenesis",
    }
    source_offsets = {"Discovery": -0.27, "External 1": -0.09, "External 2": 0.09, "Random-effects meta": 0.27}
    colors = {"Discovery": "#1f77b4", "External 1": "#2ca02c", "External 2": "#ff7f0e", "Random-effects meta": "#111111"}
    markers = {"Discovery": "o", "External 1": "s", "External 2": "^", "Random-effects meta": "D"}

    fig, ax = plt.subplots(figsize=(12, 9))
    for source in source_offsets:
        subset = data.loc[data["source"].eq(source)]
        for _, row in subset.iterrows():
            if row["feature"] not in order:
                continue
            y = order.index(row["feature"]) + source_offsets[source]
            estimate = row["estimate"]
            se = row["std_error"]
            ax.errorbar(estimate, y, xerr=1.96 * se, fmt=markers[source], color=colors[source], capsize=2.5, markersize=6 if source != "Random-effects meta" else 7, label=source if y == source_offsets[source] else None)
    ax.axvline(0, color="#616161", linewidth=1.2)
    ax.set_yticks(range(len(order)), [pretty.get(x, x) for x in order])
    ax.set_xlabel("Standardized PE-control effect (95% CI)")
    ax.set_title("Pathway effects across independent cohorts")
    handles = [plt.Line2D([0], [0], marker=markers[s], color=colors[s], linestyle="", label=s) for s in source_offsets]
    ax.legend(handles=handles, loc="lower right", frameon=True, fontsize=10)
    ax.grid(axis="y", alpha=0.15)
    save(fig, "figure2_pathway_forest.png")


def target_heatmap() -> None:
    data = pd.read_csv(TABLES / "target_gene_cohort_effects.csv")
    pivot = data.pivot(index="feature", columns="cohort", values="estimate")
    pivot = pivot.rename(columns={"GSE75010_BioBank157": "Discovery", "GSE190971": "External 1", "GSE204835_term": "External 2"})
    pivot = pivot.reindex(columns=["Discovery", "External 1", "External 2"])
    meta = pd.read_csv(TABLES / "target_gene_meta_analysis.csv").set_index("feature")["estimate"]
    pivot["Random-effects meta"] = meta
    pivot = pivot.loc[pivot["Random-effects meta"].sort_values().index]
    vmax = np.nanpercentile(np.abs(pivot.to_numpy()), 95)
    fig, ax = plt.subplots(figsize=(9, 9))
    sns.heatmap(pivot, cmap="vlag", center=0, vmin=-vmax, vmax=vmax, annot=True, fmt=".2f", linewidths=0.5, cbar_kws={"label": "Standardized effect"}, ax=ax)
    ax.set_title("Prespecified mechanism-gene effects")
    ax.set_xlabel("")
    ax.set_ylabel("")
    save(fig, "figure3_target_gene_heatmap.png")


def negative_control_plot() -> None:
    full = pd.read_csv(TABLES / "negative_control_iterations.csv.gz")
    summary = pd.read_csv(TABLES / "negative_control_summary.csv").set_index("gene_set")
    selected = [
        "HALLMARK_HYPOXIA",
        "HALLMARK_TNF_ALPHA_SIGNALING_VIA_NF_KB",
        "CURATED_ANTIANGIOGENIC_VASCULAR_STRESS",
    ]
    titles = ["Hypoxia", "TNF-alpha / NF-kB", "Antiangiogenic stress"]
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.7), sharey=True)
    for ax, gene_set, title in zip(axes, selected, titles):
        values = full.loc[(full["level_0"].eq("matched_random")) & full["gene_set"].eq(gene_set), "estimate"].dropna()
        observed = summary.loc[gene_set, "observed_estimate"]
        p = summary.loc[gene_set, "matched_random_plus_one_p"]
        sns.histplot(values, bins=35, color="#90a4ae", edgecolor="white", ax=ax)
        ax.axvline(observed, color="#c62828", linewidth=2.5)
        ax.set_title(title)
        ax.set_xlabel("")
        ax.text(0.04, 0.92, f"Observed = {observed:.2f}\nplus-one p = {p:.3f}", transform=ax.transAxes, va="top", fontsize=10, bbox=dict(facecolor="white", alpha=0.9, edgecolor="#bdbdbd"))
    axes[0].set_ylabel("Random gene sets")
    axes[1].set_xlabel("Effect from expression-matched random gene set", fontsize=13, labelpad=18)
    fig.suptitle("Matched-gene-set negative controls in the discovery cohort", fontsize=18, weight="bold", y=1.04)
    fig.subplots_adjust(bottom=0.22)
    save(fig, "figure4_negative_controls.png")


def mirna_plot() -> None:
    data = pd.read_csv(TABLES / "GSE177049_mirna_effects.csv")
    data = data.sort_values("estimate")
    fig, ax = plt.subplots(figsize=(8.5, 4.8))
    y = np.arange(len(data))
    ax.errorbar(data["estimate"], y, xerr=1.96 * data["std_error"], fmt="o", color="#6a3d9a", capsize=4)
    ax.axvline(0, color="#616161", linewidth=1.2)
    ax.set_yticks(y, data["feature"])
    ax.set_xlabel("Standardized PE-control effect (95% CI)")
    ax.set_title("Exploratory miRNA effects: GSE177049 (5 vs 5)")
    save(fig, "figure5_exploratory_mirna.png")


def volcano_plot() -> None:
    data = pd.read_csv(TABLES / "discovery_gene_differential_expression.csv.gz")
    data["minus_log10_p"] = -np.log10(data["p_value"].clip(lower=np.finfo(float).tiny))
    significant = data["fdr"] < 0.05
    fig, ax = plt.subplots(figsize=(10.5, 7))
    ax.scatter(data.loc[~significant, "estimate"], data.loc[~significant, "minus_log10_p"], s=8, color="#b0bec5", alpha=0.45, linewidths=0)
    ax.scatter(data.loc[significant, "estimate"], data.loc[significant, "minus_log10_p"], s=9, color="#1565c0", alpha=0.55, linewidths=0)
    targets = {"FLT1", "ENG", "SERPINE1", "PGF", "NOS3", "HIF1A", "TWIST1", "NFKB1", "VEGFA"}
    labelled = data["feature"].isin(targets) & data["fdr"].lt(0.05)
    for _, row in data.loc[labelled].iterrows():
        ax.scatter(row["estimate"], row["minus_log10_p"], s=35, color="#c62828", zorder=4)
        ax.annotate(row["feature"], (row["estimate"], row["minus_log10_p"]), xytext=(4, 4), textcoords="offset points", fontsize=9)
    ax.axvline(0, color="#616161", linewidth=1)
    ax.set_xlabel("Adjusted PE-control log2-expression coefficient")
    ax.set_ylabel("-log10(p-value)")
    ax.set_title("Discovery differential expression: GSE75010 BioBank 157")
    ax.text(0.02, 0.97, "Adjusted for gestational age, fetal sex, and chronic hypertension", transform=ax.transAxes, va="top", fontsize=9, bbox=dict(facecolor="white", alpha=0.85, edgecolor="none"))
    save(fig, "figure6_discovery_volcano.png")


def main() -> None:
    study_design()
    pathway_forest()
    target_heatmap()
    negative_control_plot()
    mirna_plot()
    volcano_plot()
    for path in sorted(FIGURES.glob("*.png")):
        print(f"{path.name}: {path.stat().st_size:,} bytes")


if __name__ == "__main__":
    main()
