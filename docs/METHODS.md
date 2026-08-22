# Technical methods and interpretation guide

This document explains the statistical workflow behind the retained results. It is intended for reviewers who want more detail than the main README without reading every implementation file.

## 1. Study objective

The primary question is whether placental transcriptomic programs related to hypoxia, inflammation, and angiogenic dysfunction show reproducible PE-control differences across independently collected cohorts.

The study evaluates association and cross-cohort reproducibility. It does not estimate causal effects, treatment response, diagnostic calibration, or clinical decision benefit.

## 2. Cohort-specific expression values

Let `x_ig` denote the normalized expression of gene `g` in sample `i`.

- Microarray cohorts use processed, background-corrected, normalized log-scale expression.
- RNA-seq count matrices are filtered, normalized with median-ratio size factors, and transformed as `log2(normalized_count + 0.5)`.
- Duplicate gene symbols are collapsed by their within-sample mean before analysis.
- Raw expression values are never pooled across assay platforms.

## 3. Gene and pathway standardization

Within each cohort, each gene is standardized across samples:

```text
z_ig = (x_ig - mean_g) / sd_g
```

`mean_g` and `sd_g` are the mean and sample standard deviation of gene `g` across the samples in that cohort. A value of `z_ig = 1` therefore means that sample `i` expresses gene `g` one cohort standard deviation above the cohort mean.

For pathway `P`, the sample-level pathway score is the mean standardized expression across measured member genes:

```text
score_iP = mean(z_ig for g in P)
```

At least three measurable member genes are required. The resulting pathway score is standardized once more across samples before regression, making the PE coefficient interpretable as an adjusted difference in pathway-score standard-deviation units.

## 4. Cohort models

The discovery model is:

```text
standardized_expression_or_score
    ~ intercept + PE + gestational_age + fetal_sex + chronic_hypertension
```

The coefficient on `PE` is the adjusted PE-control difference. The primary external models use `intercept + PE` because equivalent covariate fields are not consistently complete across public cohorts. Gestational-age and case-definition sensitivity analyses are retained separately where the metadata permit them.

Ordinary least squares is fitted using matrix algebra. Each result retains the estimate, classical standard error, test statistic, p-value, residual degrees of freedom, and Benjamini-Hochberg FDR.

## 5. Random-effects meta-analysis

Only independently collected cohorts enter the primary meta-analysis. For cohort estimate `y_j` with standard error `se_j`, the DerSimonian-Laird procedure estimates between-cohort variance `tau2` and applies random-effects weights:

```text
w_j = 1 / (se_j^2 + tau2)
pooled_effect = sum(w_j * y_j) / sum(w_j)
```

The pooled effect is a weighted standardized PE-control difference. It is not a similarity percentage. Heterogeneity is summarized by `tau2`, Cochran's Q, and I2.

## 6. Multiple-testing control

Benjamini-Hochberg FDR is applied within each reported testing family. FDR controls the expected proportion of false discoveries among declared discoveries under the procedure's assumptions. It does not mean that an individual result has a 5% probability of being false.

## 7. Dataset genealogy audit

GSE75010 contains 157 newly collected BioBank profiles and 173 profiles imported from seven earlier studies. The imported profiles are retained as historical internal replication but excluded from external meta-analysis. This prevents the same participants from being counted once in discovery and again as apparently independent validation.

Sample role, source study, platform, inclusion status, and provenance notes are recorded in `data/processed/cohort_manifest.csv`.

## 8. Negative controls

### 8.1 Expression/variance-matched random gene sets

For each measured pathway, 500 random gene sets are generated with approximately the same number of measurable genes. Candidate genes are matched in ten-by-ten bins of cohort mean expression and expression standard deviation. Real pathway members and duplicate selections are excluded.

Each random set is scored and modeled with the same pipeline as the real pathway. The two-sided plus-one empirical p-value is:

```text
p = (1 + number of |random effects| >= |observed effect|) / (B + 1)
```

This control asks whether the named pathway is unusually strong relative to arbitrary gene sets with comparable basic expression properties.

### 8.2 Outcome-label permutation

The observed pathway scores, sample covariates, and pathway membership are held fixed. PE-control labels are permuted 500 times, and the PE coefficient is refitted after every permutation.

Patient scores do not change. What changes is which scores are assigned to the permuted PE and control groups. This control asks whether association with the observed diagnostic labels is stronger than expected under random label assignments.

With 500 iterations, zero null effects as extreme as the observed effect yields `(0 + 1) / (500 + 1) = 0.001996`, not zero.

The matched-set and label-permutation controls answer different questions. A pathway can strongly separate the observed labels without being uniquely stronger than all comparable random gene sets.

## 9. Leave-one-cohort-out analysis

The independent meta-analysis is repeated after omitting each cohort in turn:

1. omit GSE75010 BioBank and combine the two external cohorts;
2. omit GSE190971 and combine discovery with GSE204835;
3. omit GSE204835 and combine discovery with GSE190971.

Stable direction and similar magnitude across these analyses indicate that no single cohort entirely determines the pooled result. With only two cohorts remaining in each iteration, this is a sensitivity analysis rather than proof of universal generalizability.

## 10. Exploratory paired mRNA/miRNA analysis

GSE177049 provides mature-miRNA and mRNA measurements from the same ten placental samples. This pairing permits:

- standardized PE-control comparisons for selected mature miRNAs; and
- within-sample Spearman correlations between selected miRNAs and target mRNAs.

Most mRNA cohorts do not directly measure mature miRNA, and host-transcript abundance is not a reliable substitute for mature-miRNA abundance. The paired design supports association analysis but not causality. The cohort is small, measurements are cross-sectional, and correlations may reflect disease severity, cell composition, shared regulation, or chance. No examined miRNA or prespecified miRNA-target correlation passed FDR 0.05.

## 11. Reproducibility artifacts

The repository retains:

- official source URLs and expected source hashes;
- sample-level cohort genealogy and inclusion decisions;
- frozen analysis settings and prespecified gene sets;
- cohort-level estimates and pathway scores;
- meta-analysis and leave-one-cohort-out tables;
- all matched-random and label-permutation iterations;
- exploratory miRNA effects and correlations;
- figures and a generated PDF report.

Raw source files and processed expression matrices are rebuilt locally and excluded from Git because they are downloadable and comparatively large.
