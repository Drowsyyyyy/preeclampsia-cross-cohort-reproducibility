"""Build the portfolio-ready PDF technical report from retained artifacts."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from reportlab.lib import colors
from reportlab.lib.colors import HexColor
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import inch
from reportlab.platypus import (
    HRFlowable,
    Image,
    PageBreak,
    Paragraph,
    SimpleDocTemplate,
    Spacer,
    Table,
    TableStyle,
)


ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "results" / "tables"
FIGURES = ROOT / "results" / "figures"
OUTPUT = ROOT / "output" / "pdf" / "preeclampsia_cross_cohort_report.pdf"

NAVY = HexColor("#17324d")
BLUE = HexColor("#2b6f9f")
LIGHT_BLUE = HexColor("#e9f2f8")
GREEN = HexColor("#2d7f5e")
LIGHT_GREEN = HexColor("#eaf4ef")
ORANGE = HexColor("#9a5b13")
LIGHT_ORANGE = HexColor("#fff3df")
PURPLE = HexColor("#68469b")
LIGHT_GREY = HexColor("#f4f6f7")
DARK_GREY = HexColor("#333b42")


def clean(text: str) -> str:
    """Keep PDF copy text portable and use ASCII hyphens."""
    return (
        str(text)
        .replace("\u2013", "-")
        .replace("\u2014", "-")
        .replace("\u2212", "-")
        .replace("\u00a0", " ")
    )


def fmt_p(value: float) -> str:
    if value < 0.001:
        return f"{value:.2e}"
    return f"{value:.3f}"


def add_page_number(canvas, doc) -> None:
    canvas.saveState()
    canvas.setStrokeColor(HexColor("#d8dde1"))
    canvas.line(0.65 * inch, 0.55 * inch, 7.85 * inch, 0.55 * inch)
    canvas.setFont("Helvetica", 8)
    canvas.setFillColor(HexColor("#6b747c"))
    canvas.drawString(0.68 * inch, 0.35 * inch, "Preeclampsia cross-cohort reproducibility study")
    canvas.drawRightString(7.82 * inch, 0.35 * inch, f"Page {doc.page}")
    canvas.restoreState()


def scaled_image(path: Path, max_width: float, max_height: float) -> Image:
    from PIL import Image as PILImage

    with PILImage.open(path) as picture:
        width, height = picture.size
    scale = min(max_width / width, max_height / height)
    return Image(str(path), width=width * scale, height=height * scale)


def make_styles():
    base = getSampleStyleSheet()
    return {
        "title": ParagraphStyle(
            "Title",
            parent=base["Title"],
            fontName="Helvetica-Bold",
            fontSize=25,
            leading=30,
            textColor=NAVY,
            alignment=TA_LEFT,
            spaceAfter=14,
        ),
        "subtitle": ParagraphStyle(
            "Subtitle",
            parent=base["Normal"],
            fontName="Helvetica",
            fontSize=12,
            leading=18,
            textColor=DARK_GREY,
            spaceAfter=12,
        ),
        "h1": ParagraphStyle(
            "H1",
            parent=base["Heading1"],
            fontName="Helvetica-Bold",
            fontSize=17,
            leading=21,
            textColor=NAVY,
            spaceBefore=5,
            spaceAfter=10,
        ),
        "h2": ParagraphStyle(
            "H2",
            parent=base["Heading2"],
            fontName="Helvetica-Bold",
            fontSize=12.5,
            leading=16,
            textColor=BLUE,
            spaceBefore=8,
            spaceAfter=5,
        ),
        "body": ParagraphStyle(
            "Body",
            parent=base["BodyText"],
            fontName="Helvetica",
            fontSize=9.4,
            leading=13.2,
            textColor=DARK_GREY,
            spaceAfter=7,
        ),
        "small": ParagraphStyle(
            "Small",
            parent=base["BodyText"],
            fontName="Helvetica",
            fontSize=7.8,
            leading=10.5,
            textColor=HexColor("#505960"),
            spaceAfter=4,
        ),
        "caption": ParagraphStyle(
            "Caption",
            parent=base["BodyText"],
            fontName="Helvetica-Oblique",
            fontSize=7.8,
            leading=10.5,
            textColor=HexColor("#505960"),
            spaceBefore=5,
            spaceAfter=8,
        ),
        "callout": ParagraphStyle(
            "Callout",
            parent=base["BodyText"],
            fontName="Helvetica",
            fontSize=10.2,
            leading=14.5,
            textColor=NAVY,
            borderColor=BLUE,
            borderWidth=1,
            borderPadding=10,
            backColor=LIGHT_BLUE,
            spaceAfter=10,
        ),
        "bullet": ParagraphStyle(
            "Bullet",
            parent=base["BodyText"],
            fontName="Helvetica",
            fontSize=9.2,
            leading=12.8,
            leftIndent=13,
            firstLineIndent=-8,
            bulletIndent=3,
            textColor=DARK_GREY,
            spaceAfter=5,
        ),
        "center_small": ParagraphStyle(
            "CenterSmall",
            parent=base["BodyText"],
            fontName="Helvetica",
            fontSize=8.2,
            leading=11,
            alignment=TA_CENTER,
            textColor=HexColor("#505960"),
        ),
    }


def P(text: str, style) -> Paragraph:
    return Paragraph(clean(text), style)


def data_table(rows, widths, styles, header=True, font_size=7.7) -> Table:
    converted = []
    for row_index, row in enumerate(rows):
        row_style = styles["small"]
        if header and row_index == 0:
            row_style = ParagraphStyle(
                "TableHeader",
                parent=styles["small"],
                fontName="Helvetica-Bold",
                textColor=colors.white,
                alignment=TA_CENTER,
                fontSize=font_size,
            )
        converted.append([P(str(cell), row_style) for cell in row])
    table = Table(converted, colWidths=widths, repeatRows=1 if header else 0, hAlign="LEFT")
    commands = [
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("LEFTPADDING", (0, 0), (-1, -1), 5),
        ("RIGHTPADDING", (0, 0), (-1, -1), 5),
        ("TOPPADDING", (0, 0), (-1, -1), 5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
        ("GRID", (0, 0), (-1, -1), 0.35, HexColor("#cfd6db")),
    ]
    if header:
        commands.extend(
            [
                ("BACKGROUND", (0, 0), (-1, 0), NAVY),
                ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, LIGHT_GREY]),
            ]
        )
    table.setStyle(TableStyle(commands))
    return table


def main() -> None:
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    styles = make_styles()

    summary = json.loads((ROOT / "results" / "analysis_summary.json").read_text(encoding="utf-8"))
    pathway_meta = pd.read_csv(TABLES / "independent_pathway_meta_analysis.csv").set_index("feature")
    gene_meta = pd.read_csv(TABLES / "target_gene_meta_analysis.csv").set_index("feature")
    composition = pd.read_csv(TABLES / "composition_adjusted_pathways.csv").set_index("feature")
    negative = pd.read_csv(TABLES / "negative_control_summary.csv").set_index("gene_set")
    mirna = pd.read_csv(TABLES / "GSE177049_mirna_effects.csv").set_index("feature")

    document = SimpleDocTemplate(
        str(OUTPUT),
        pagesize=letter,
        rightMargin=0.65 * inch,
        leftMargin=0.65 * inch,
        topMargin=0.65 * inch,
        bottomMargin=0.72 * inch,
        title="Cross-cohort evaluation of preeclampsia transcriptomes",
        author="Jungwon Kim",
        subject="Reproducible public-data transcriptomics study",
    )
    story = []

    # Cover
    story.append(Spacer(1, 0.45 * inch))
    story.append(P("Cross-cohort evaluation of hypoxia, inflammation, and angiogenic dysfunction in human preeclampsia transcriptomes", styles["title"]))
    story.append(HRFlowable(width="100%", thickness=2.2, color=BLUE, spaceBefore=4, spaceAfter=16))
    story.append(P("A reproducibility-first public-data pilot", styles["subtitle"]))
    story.append(
        P(
            "<b>Objective.</b> Test whether placental transcriptomic patterns across independently collected human preeclampsia cohorts are compatible with prior hypoxia, inflammatory, and angiogenic mechanisms - while separating genuine external replication from samples already embedded in an integrated dataset.",
            styles["callout"],
        )
    )
    cover_rows = [
        ["Project status", "Completed computational study using public, de-identified data"],
        ["Analysis date", "2026-08-21"],
        ["Discovery", "GSE75010 new BioBank: 80 PE, 77 non-PE"],
        ["Independent replication", "GSE190971: 7 PE, 6 controls; GSE204835 term: 12 PE, 12 controls"],
        ["Exploratory miRNA", "GSE177049: 5 early-onset PE, 5 preterm controls"],
        ["Primary framework", "Within-cohort modeling, standardized effects, random-effects meta-analysis"],
    ]
    cover = data_table(cover_rows, [1.55 * inch, 5.15 * inch], styles, header=False, font_size=8)
    cover.setStyle(TableStyle([("BACKGROUND", (0, 0), (0, -1), LIGHT_BLUE), ("FONTNAME", (0, 0), (0, -1), "Helvetica-Bold")]))
    story.append(cover)
    story.append(Spacer(1, 0.3 * inch))
    story.append(P("Prepared by Jungwon Kim", styles["h2"]))
    story.append(P("This document reports a current independent analysis. It does not represent a historical collaboration, faculty-supervised project, peer-reviewed publication, clinical validation, or causal study.", styles["small"]))
    story.append(PageBreak())

    # Executive summary
    story.append(P("Executive summary", styles["h1"]))
    story.append(
        P(
            "The clearest reproducible signals were <b>hypoxia</b>, <b>antiangiogenic/vascular stress</b>, and <b>TNF-alpha/NF-kB signaling</b>. Their standardized PE-control effects were positive in the discovery cohort and both independent replication cohorts, with random-effects meta-analysis FDR values below 0.05. FLT1, ENG, and SERPINE1 were among the most consistent prespecified genes.",
            styles["callout"],
        )
    )
    bullets = [
        f"Discovery differential expression identified {summary['discovery_significant_genes_fdr_0_05']:,} genes at FDR 0.05; this large number is descriptive, not a claim that every gene is mechanistically important.",
        f"Independent meta-analysis identified {summary['independent_meta_significant_pathways_fdr_0_05']} of 10 measured pathways and {summary['independent_meta_significant_genes_fdr_0_05']:,} genes at FDR 0.05.",
        "Hypoxia and antiangiogenic stress each had a pooled standardized effect near 0.96 with no estimated pathway-level heterogeneity in the three primary cohorts.",
        "PGF, NOS3, and broader proangiogenic/NO signaling were heterogeneous. The data do not support a universal direction across cohorts.",
        "The exploratory miRNA cohort was very small (5 versus 5); no miRNA or miRNA-target correlation passed FDR 0.05.",
        "Matched random gene sets and label permutations were used as negative controls. Hypoxia was unusually strong relative to expression-matched random sets; TNF/NF-kB and antiangiogenic stress separated labels but were not uniquely extreme under the matched-set comparison.",
    ]
    for item in bullets:
        story.append(P(f"&#8226; {item}", styles["bullet"]))
    story.append(P("Bottom line", styles["h2"]))
    story.append(P("The study supports a reproducible placental stress signature, but it does not show that a specific miRNA causes preeclampsia, that bulk mRNA identifies the responsible cell type, or that these signals have diagnostic or treatment value.", styles["body"]))
    story.append(PageBreak())

    # Design and provenance
    story.append(P("1. Evidence hierarchy and data provenance", styles["h1"]))
    story.append(scaled_image(FIGURES / "figure1_study_design.png", 7.1 * inch, 4.2 * inch))
    story.append(P("Figure 1. Evidence hierarchy. The seven historical studies are already incorporated into GSE75010, so they are used only for internal historical replication and are never counted as external validation.", styles["caption"]))
    story.append(P("The central design decision was to audit dataset genealogy before analysis. GSE75010 contains 330 profiles: 157 newly collected RCWIH BioBank samples and 173 selected PE/control profiles imported from seven earlier studies. Calling those seven studies external validation would reuse the same profiles twice and overstate reproducibility.", styles["body"]))
    story.append(P("The independent evidence therefore consists of the new GSE75010 BioBank samples plus two separately collected GEO cohorts. Every cohort was processed on its own scale. Only standardized within-cohort effects were combined.", styles["body"]))
    story.append(P("Recorded discrepancy", styles["h2"]))
    story.append(P("The GSE204835 series summary reports 11 term and 5 preterm de novo PE samples, while current sample-level SOFT fields identify 12 term and 4 preterm samples. The primary analysis follows the sample-level fields and records this discrepancy in data/processed/data_manifest.json.", styles["body"]))
    story.append(PageBreak())

    # Methods
    story.append(P("2. Methods", styles["h1"]))
    methods = [
        ("Prespecification", "The cohort roles, contrasts, target genes, pathway definitions, covariates, FDR threshold, meta-analysis method, negative controls, and interpretation limits were written to config/analysis_plan.yaml before outcome analysis."),
        ("Discovery model", "For GSE75010 BioBank samples, gene-level log2 expression was modeled as PE status plus gestational age, fetal sex, and chronic hypertension. Duplicate gene symbols had already been collapsed during preparation."),
        ("RNA-seq cohorts", "Counts were normalized within each cohort using median-ratio size factors and log2(normalized count + 0.5). GSE190971 compared all placental PE with controls. The GSE204835 primary contrast used term de novo PE and term controls."),
        ("Pathway score", "For each measured gene set, genes were standardized within cohort and averaged per sample. Group effects were converted to standardized coefficients with standard errors."),
        ("Meta-analysis", "DerSimonian-Laird random-effects models combined within-cohort standardized effects. No raw expression values were pooled across assay platforms."),
        ("Multiplicity", "Benjamini-Hochberg FDR controlled the reported gene- and pathway-level multiple testing families."),
        ("Robustness", "The study retained leave-one-cohort-out meta-analysis, discovery models adjusted for placental cell-marker proxy scores, 500 expression/variance-matched random gene sets, and 500 PE-label permutations per measured pathway."),
        ("Exploratory miRNA", "GSE177049 mature-miRNA effects and Spearman correlations with selected mRNA targets were descriptive because the cohort contains only 10 samples."),
    ]
    method_rows = [["Component", "Implementation"]] + [[name, text] for name, text in methods]
    story.append(data_table(method_rows, [1.45 * inch, 5.25 * inch], styles))
    story.append(Spacer(1, 8))
    story.append(P("Effect direction", styles["h2"]))
    story.append(P("Positive standardized effects indicate higher expression or pathway scores in PE than in controls. Confidence intervals that cross zero are not interpreted as statistically resolved. FDR is emphasized over nominal p-values.", styles["body"]))
    story.append(PageBreak())

    # Pathways
    story.append(P("3. Independent pathway replication", styles["h1"]))
    story.append(scaled_image(FIGURES / "figure2_pathway_forest.png", 7.1 * inch, 4.7 * inch))
    story.append(P("Figure 2. Standardized pathway effects and 95% confidence intervals in the three independently collected cohorts. Black diamonds are random-effects pooled estimates. Some curated sets were unavailable when fewer than three genes were measured.", styles["caption"]))
    selected_paths = [
        ("CURATED_ANTIANGIOGENIC_VASCULAR_STRESS", "Antiangiogenic/vascular stress"),
        ("HALLMARK_HYPOXIA", "Hypoxia"),
        ("HALLMARK_TNF_ALPHA_SIGNALING_VIA_NF_KB", "TNF-alpha/NF-kB"),
        ("CURATED_HIF1_TWIST1_MIR214_HOST_AXIS", "HIF1/TWIST1/miR-214 host"),
        ("HALLMARK_INFLAMMATORY_RESPONSE", "Inflammatory response"),
        ("CURATED_PROANGIOGENIC_SUPPORT", "Proangiogenic support"),
        ("CURATED_ENDOTHELIAL_NO_SIGNALING", "Endothelial NO signaling"),
    ]
    pathway_rows = [["Pathway", "Effect (95% CI)", "FDR", "I2", "k"]]
    for key, label in selected_paths:
        row = pathway_meta.loc[key]
        pathway_rows.append([label, f"{row.estimate:.2f} ({row.ci_low:.2f}, {row.ci_high:.2f})", fmt_p(row.fdr), f"{row.i2:.0f}%", f"{int(row.k)}"])
    story.append(data_table(pathway_rows, [2.55 * inch, 1.55 * inch, 0.85 * inch, 0.7 * inch, 0.45 * inch], styles))
    story.append(PageBreak())

    # Target genes
    story.append(P("4. Prespecified mechanism genes", styles["h1"]))
    story.append(scaled_image(FIGURES / "figure3_target_gene_heatmap.png", 4.6 * inch, 4.9 * inch))
    story.append(P("Figure 3. Cohort-specific and random-effects standardized effects for prespecified genes. Blank cells indicate that the gene was not reliably available in that cohort. Color shows direction and magnitude, not statistical significance.", styles["caption"]))
    important_genes = ["FLT1", "ENG", "SERPINE1", "NFKB1", "TWIST1", "HIF1A", "NOS3", "PGF"]
    gene_rows = [["Gene", "Meta effect (95% CI)", "FDR", "I2", "Interpretive note"]]
    notes = {
        "FLT1": "Strong positive; gene-level signal does not isolate soluble sFlt-1.",
        "ENG": "Consistent positive vascular-stress marker.",
        "SERPINE1": "Positive with no estimated heterogeneity.",
        "NFKB1": "Modest positive transcript-level support.",
        "TWIST1": "Nominal negative, not FDR-significant.",
        "HIF1A": "Not significant; mRNA is not HIF-1alpha protein activity.",
        "NOS3": "Not significant; measured in only two cohorts.",
        "PGF": "Highly heterogeneous, including opposite direction in external 2.",
    }
    for gene in important_genes:
        row = gene_meta.loc[gene]
        gene_rows.append([gene, f"{row.estimate:.2f} ({row.ci_low:.2f}, {row.ci_high:.2f})", fmt_p(row.fdr), f"{row.i2:.0f}%", notes[gene]])
    story.append(data_table(gene_rows, [0.8 * inch, 1.35 * inch, 0.65 * inch, 0.55 * inch, 3.35 * inch], styles, font_size=6.8))
    story.append(PageBreak())

    # Robustness
    story.append(P("5. Robustness and negative controls", styles["h1"]))
    story.append(scaled_image(FIGURES / "figure4_negative_controls.png", 7.1 * inch, 3.2 * inch))
    story.append(P("Figure 4. Discovery effects (red) compared with 500 random gene sets matched approximately on size, expression, and variability. The plus-one proportion is (extreme random sets + 1)/(500 + 1).", styles["caption"]))
    neg_rows = [["Pathway", "Observed", "Matched-set plus-one p", "Label-permutation plus-one p"]]
    for key, label in [
        ("HALLMARK_HYPOXIA", "Hypoxia"),
        ("HALLMARK_TNF_ALPHA_SIGNALING_VIA_NF_KB", "TNF-alpha/NF-kB"),
        ("CURATED_ANTIANGIOGENIC_VASCULAR_STRESS", "Antiangiogenic stress"),
        ("HALLMARK_INFLAMMATORY_RESPONSE", "Inflammatory response"),
    ]:
        row = negative.loc[key]
        neg_rows.append([label, f"{row.observed_estimate:.2f}", fmt_p(row.matched_random_plus_one_p), fmt_p(row.label_permutation_plus_one_p)])
    story.append(data_table(neg_rows, [2.3 * inch, 1.0 * inch, 1.7 * inch, 1.7 * inch], styles))
    story.append(P("Interpretation", styles["h2"]))
    story.append(P("Hypoxia was more extreme than almost all matched random sets. The antiangiogenic and TNF/NF-kB scores were strongly associated with the observed PE labels, but their magnitudes were not uniquely extreme relative to all matched random sets at the 0.05 threshold. This limits pathway-specific claims even when label permutations are significant.", styles["body"]))
    story.append(P("Cell-composition sensitivity", styles["h2"]))
    comp_items = []
    for key, label in [
        ("CURATED_ANTIANGIOGENIC_VASCULAR_STRESS", "Antiangiogenic stress"),
        ("HALLMARK_HYPOXIA", "Hypoxia"),
        ("HALLMARK_TNF_ALPHA_SIGNALING_VIA_NF_KB", "TNF-alpha/NF-kB"),
    ]:
        row = composition.loc[key]
        comp_items.append(f"{label}: adjusted effect {row.estimate:.2f}, FDR {fmt_p(row.fdr)}")
    story.append(P("; ".join(comp_items) + ". These marker-score adjustments are sensitivity analyses, not formal cell-type deconvolution.", styles["body"]))
    story.append(PageBreak())

    # miRNA
    story.append(P("6. Exploratory paired mRNA/miRNA analysis", styles["h1"]))
    story.append(scaled_image(FIGURES / "figure5_exploratory_mirna.png", 6.8 * inch, 4.2 * inch))
    story.append(P("Figure 5. Exploratory standardized miRNA effects in GSE177049. None passed FDR 0.05.", styles["caption"]))
    mirna_rows = [["miRNA", "Effect (95% CI)", "Nominal p", "FDR"]]
    for feature, row in mirna.sort_values("estimate").iterrows():
        low = row.estimate - 1.96 * row.std_error
        high = row.estimate + 1.96 * row.std_error
        mirna_rows.append([feature, f"{row.estimate:.2f} ({low:.2f}, {high:.2f})", fmt_p(row.p_value), fmt_p(row.fdr)])
    story.append(data_table(mirna_rows, [1.7 * inch, 2.1 * inch, 1.2 * inch, 1.0 * inch], styles))
    story.append(Spacer(1, 8))
    story.append(P("miR-155-5p had a nominal negative effect (p = 0.027) but did not survive FDR correction. None of the prespecified miRNA-target Spearman correlations passed FDR 0.05. These observations cannot be used as mechanistic confirmation.", styles["body"]))
    story.append(PageBreak())

    # Discovery context
    story.append(P("7. Discovery-wide context", styles["h1"]))
    story.append(scaled_image(FIGURES / "figure6_discovery_volcano.png", 6.8 * inch, 5.0 * inch))
    story.append(P("Figure 6. Discovery gene-level results. Prespecified genes are labelled only when they passed discovery FDR 0.05. The y-axis uses nominal p-values for display; formal interpretation uses FDR.", styles["caption"]))
    story.append(P("The discovery cohort contains broad transcriptomic separation between PE and non-PE placentas. FLT1 is the strongest prespecified signal, followed by ENG and SERPINE1. The volcano plot is not evidence that every statistically significant gene is a biomarker or causal mediator; gestational age, treatment, placental composition, disease severity, and other factors may influence bulk expression.", styles["body"]))
    story.append(PageBreak())

    # Limitations and conclusion
    story.append(P("8. Interpretation, limitations, and responsible use", styles["h1"]))
    limitations = [
        "All data are observational, cross-sectional placental profiles collected at delivery. Directionality and causality cannot be established.",
        "GSE190971 is very small, and the primary GSE204835 comparison contains 24 samples. Confidence intervals are therefore wide for many pathways.",
        "Bulk tissue mixes trophoblast, endothelial, immune, stromal, and other cell populations. Marker-score adjustment cannot replace single-cell or spatial validation.",
        "HIF1A transcript abundance is not a direct measurement of HIF-1alpha stabilization, DNA binding, or transcriptional activity.",
        "Gene-level FLT1 cannot distinguish all transcript isoforms and does not directly quantify circulating sFlt-1 protein.",
        "Mature miR-31, miR-155, and miR-214 are absent from most mRNA cohorts. The one paired miRNA cohort is underpowered.",
        "Platform, gestational-age distribution, disease definition, tissue processing, and ancestry differences limit cross-cohort equivalence.",
        "Discrimination, calibration, clinical decision benefit, treatment response, and deployment-laboratory robustness were not assessed.",
    ]
    for item in limitations:
        story.append(P(f"&#8226; {item}", styles["bullet"]))
    story.append(P("Conclusion", styles["h2"]))
    story.append(
        P(
            "Across three independently collected placental cohorts, preeclampsia was consistently associated with hypoxia, antiangiogenic/vascular stress, and TNF-alpha/NF-kB transcriptomic programs. FLT1, ENG, and SERPINE1 provided strong gene-level support. Other hypothesized elements - especially PGF/NOS3 direction and miRNA regulation - were heterogeneous or unresolved. The contribution is therefore methodological and evidential: a transparent workflow that distinguishes robust replication from attractive but insufficiently supported mechanistic claims.",
            styles["callout"],
        )
    )
    story.append(PageBreak())

    # Reproducibility and references
    story.append(P("9. Reproducibility record", styles["h1"]))
    reproducibility = [
        ("Frozen plan", "config/analysis_plan.yaml"),
        ("Source URLs", "config/source_urls.yaml"),
        ("Download", "scripts/00_download_data.py"),
        ("Preparation", "scripts/01_prepare_data.py"),
        ("Analysis", "scripts/02_run_analysis.py"),
        ("Figures", "scripts/03_make_figures.py"),
        ("Report", "scripts/04_build_report.py"),
        ("Quality checks", "scripts/99_quality_checks.py"),
        ("Data genealogy", "data/processed/cohort_manifest.csv"),
        ("Hashes and discrepancies", "data/processed/data_manifest.json"),
        ("Fold-free retained estimates", "results/tables/ and results/cohort_level_results/"),
    ]
    story.append(data_table([["Record", "Location"]] + [[a, b] for a, b in reproducibility], [2.0 * inch, 4.7 * inch], styles))
    story.append(Spacer(1, 10))
    story.append(P("The pipeline stores source-file SHA-256 hashes, every sample's inclusion role, cohort-level estimates, meta-analysis inputs, negative-control iterations, final tables, and figures. Running scripts/99_quality_checks.py verifies hashes, expected sample counts, non-overlap of embedded and discovery profiles, matrix finiteness, independent-cohort restrictions, control iteration counts, and retained outputs.", styles["body"]))
    story.append(P("Primary data resources", styles["h2"]))
    resources = [
        ("GSE75010", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75010"),
        ("GSE190971", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE190971"),
        ("GSE204835", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE204835"),
        ("GSE177049", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE177049"),
        ("NCBI Gene", "https://www.ncbi.nlm.nih.gov/gene/"),
        ("Enrichr libraries", "https://maayanlab.cloud/Enrichr/"),
    ]
    for name, url in resources:
        story.append(P(f"&#8226; <b>{name}</b>: <link href='{url}' color='#2b6f9f'>{url}</link>", styles["small"]))
    story.append(P("Ethics and scope", styles["h2"]))
    story.append(P("Only public, de-identified secondary data were used. No attempt was made to re-identify participants. This analysis is not medical advice and is not intended for patient-level prediction or clinical decision making.", styles["body"]))

    document.build(story, onFirstPage=add_page_number, onLaterPages=add_page_number)
    print(f"created {OUTPUT} ({OUTPUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
