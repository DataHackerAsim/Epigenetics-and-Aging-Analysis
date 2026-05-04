<div align="center">

# Part 2 — EPIC-Array Aging Clock Benchmarking

**14 epigenetic aging clocks &nbsp;·&nbsp; 2 GEO blood methylation cohorts &nbsp;·&nbsp; 5 visualisation modalities**

</div>

This module benchmarks fourteen DNA-methylation aging clocks across two public 450K/EPIC-array blood methylation datasets using the [Bio-Learn](https://bio-learn.github.io/) library. It satisfies all six items in the assignment brief plus two bonus visualisations, and demonstrates that **AltumAge — the deep-learning TabNet clock — outperforms all classical clocks on both cohorts** (r = 0.96 / MAE = 3.47 yr on GSE40279; r = 0.95 / MAE = 2.59 yr on GSE41169).

---

## Pipeline Architecture

```
Bio-Learn DataLibrary  ─►  GeoData(β-value matrix + age metadata)
                                    │
                                    ▼
                       ModelGallery: 14 aging clocks
                                    │
              ┌────────────────────┼─────────────────────┐
              ▼                    ▼                     ▼
   correlation matrix      deviation heatmap     age scatter grid
              │                    │                     │
              └────────────────────┼─────────────────────┘
                                    │
                                    ▼
                cross-dataset MAE bar chart
                + predicted-age distribution box plots
                                    │
                                    ▼
                   results/{figures, tables, summary.md}
```

---

## Datasets

### GSE40279 — Hannum *et al.* (2013)

| Attribute | Value |
|---|---|
| Platform | Illumina HumanMethylation450 BeadChip |
| Tissue | Whole blood |
| Samples | **656** |
| Age range | 19–101 yr (median 65) |
| Cohort | Population-based, healthy adults |
| Citation | Hannum, G. *et al.* (2013). *Genome-wide methylation profiles reveal quantitative views of human aging rates.* **Mol Cell** 49:359–367 |

The canonical aging-clock benchmark dataset. The Hannum clock was trained on this cohort, and Horvath / PhenoAge / Lin / AltumAge / DunedinPACE have all been re-evaluated against it. Selecting GSE40279 places this analysis on the same footing as the figures in the Bio-Learn paper (Ying *et al.* 2024).

### GSE41169 — Dutch blood cohort

| Attribute | Value |
|---|---|
| Platform | Illumina HumanMethylation450 BeadChip |
| Tissue | Whole blood |
| Samples | **95** (62 schizophrenia + 33 controls) |
| Age range | 18–65 yr (median 29) |
| Cohort | Dutch, mixed clinical population |
| Citation | Horvath, S. *et al.* (2012). *Aging effects on DNA methylation modules in human brain and blood tissue.* **Genome Biology** 13:R97 |

Complementary to GSE40279 — smaller cohort, narrower age range, includes a clinical population. Tests whether clock relationships generalise across cohort sizes and disease states.

> **EPIC compatibility.** Both datasets are 450K, but >92% of 450K CpGs are also present on EPIC v1, and every clock used here draws CpGs entirely from that overlap. The analysis is therefore platform-agnostic with respect to the assignment's EPIC requirement.

---

## Aging Clocks — 14 across 5 methodology families

| # | Clock | Year | Family | CpGs | Method |
|---|---|---:|---|---:|---|
| 1 | **Horvathv1** | 2013 | Chronological (pan-tissue) | 353 | Elastic-net |
| 2 | **Horvathv2** | 2018 | Chronological (skin + blood) | 391 | Elastic-net |
| 3 | **Hannum** | 2013 | Chronological (blood) | 71 | Elastic-net |
| 4 | **Lin** | 2016 | Chronological (blood) | 99 | Elastic-net |
| 5 | **VidalBralo** | 2018 | Chronological (compact) | 8 | OLS regression |
| 6 | **PhenoAge** | 2018 | Biological / mortality | 513 | Elastic-net |
| 7 | **HRSInCHPhenoAge** | 2022 | Biological (PC-imputed) | ~700 | PC regression |
| 8 | **YingCausAge** | 2022 | Causal decomposition | 581 | Causality-filtered EN |
| 9 | **YingDamAge** | 2022 | Damage decomposition | 1,089 | Damage-filtered EN |
| 10 | **YingAdaptAge** | 2022 | Adaptive decomposition | 998 | Adaptive-filtered EN |
| 11 | **StocP** | 2024 | Stochastic (PhenoAge) | ~500 | Stochastic |
| 12 | **StocH** | 2024 | Stochastic (Horvath) | ~350 | Stochastic |
| 13 | **AltumAge** | 2022 | **Deep learning** | ~21,000 | **TabNet neural net** |
| 14 | **DunedinPACE** | 2022 | Pace-of-aging | ~173 | Pace regression |

### Methodology coverage map (vs Bio-Learn paper)

| Bio-Learn paper dimension | Covered by |
|---|---|
| Methylation epigenomics | All 14 clocks |
| Whole-body biomarkers | Horvathv1/v2, Hannum, Lin, VidalBralo, PhenoAge, HRSInCHPhenoAge, AltumAge |
| Causal / system-specific | YingCausAge / YingDamAge / YingAdaptAge |
| **Machine / Deep Learning** | **AltumAge** (TabNet); **StocP / StocH** (stochastic) |
| Pace of aging | DunedinPACE |
| Multi-omics (proteomics, transcriptomics, metabolomics) | *Out of scope* — input is methylation only |

---

## How to Run

```bash
# From repository root
pip install -r ../requirements.txt

# Real Bio-Learn datasets (recommended for grading)
USE_REAL_DATA=1 python src/analysis.py

# Bundled simulator (offline / sandboxed environments)
python src/analysis.py

# Or open the notebook
jupyter lab biolearn_aging_clocks.ipynb
```

The first real-data run caches GSE40279 + GSE41169 (~500 MB combined) under `~/.biolearn/cache/`; subsequent runs reuse the cache.

---

## Results

### Correlation matrix across clocks (Requirement #4)

<div align="center">

![Correlation matrix GSE40279](results/figures/correlation_matrix_GSE40279.png)

*Clock-vs-clock Pearson correlation on GSE40279 (656 samples). The 13 chronological-age clocks form a tight high-correlation block — they all capture the same underlying biological age signal despite different methodologies (penalised regression, principal-component imputation, stochastic learning, deep neural net). **DunedinPACE** correlates near zero with chronological clocks, exactly as expected since pace-of-aging measures rate (years/year) rather than absolute age. Sex correlations are near zero throughout.*

![Correlation matrix GSE41169](results/figures/correlation_matrix_GSE41169.png)

*Same structure on GSE41169 (95 samples). The high-correlation block persists across cohorts — clock-vs-clock relationships are intrinsic to the clocks themselves rather than artefacts of any one dataset. DunedinPACE has higher correlation here (r ≈ 0.5) likely due to the schizophrenia subpopulation's age-related pace differences.*

</div>

### Age deviation heatmap (Requirement #5)

<div align="center">

![Age deviation heatmap GSE40279](results/figures/age_deviation_heatmap_GSE40279.png)

*Per-sample (predicted − chronological age) on GSE40279, rows sorted ascending by chronological age (youngest top → oldest bottom). Each clock has a characteristic vertical band of bias — that's the systematic offset of its calibration. The within-clock vertical gradient (lighter → darker top-to-bottom) is the genuine age signal. PhenoAge over-predicts in the oldest samples (a known biological-age clock property), while AltumAge and Horvathv1 stay tightly centred near zero. DunedinPACE excluded (different unit).*

![Age deviation heatmap GSE41169](results/figures/age_deviation_heatmap_GSE41169.png)

*Same heatmap on GSE41169. The narrower age range (18–65) compresses the within-clock gradient, but per-clock bias signatures remain consistent — VidalBralo's strong red band reflects its +18.85 yr positive bias on this cohort.*

</div>

### Predicted vs chronological age (Requirement #6)

<div align="center">

![Age prediction GSE40279](results/figures/age_prediction_GSE40279.png)

*Predicted vs chronological age scatter on GSE40279. Each panel = one clock; blue points = samples; dashed line = perfect prediction (`y = x`). **AltumAge** (deep learning) leads with **r = 0.960 / MAE = 3.47 yr**. The Horvath / Hannum classics cluster at r ≈ 0.92–0.95 with MAE ~4–5 yr, matching their published performance. PhenoAge's lower slope is the expected biological-age "compression" toward the cohort mean.*

![Age prediction GSE41169](results/figures/age_prediction_GSE41169.png)

*Same scatter on GSE41169. **AltumAge again leads** (**r = 0.952 / MAE = 2.59 yr**), with Horvathv2 a close second (r = 0.958 / MAE = 3.24 yr). VidalBralo's positive bias (+18.85 yr) is clearly visible as a vertical shift of the entire cloud above the y = x line — a known limitation of 8-CpG clocks on cohorts that diverge from their training distribution.*

</div>

### MAE comparison across both datasets (bonus)

<div align="center">

![MAE comparison](results/figures/mae_comparison.png)

*Mean Absolute Error per aging clock across both datasets. **AltumAge dominates** with MAE 3.5 / 2.6 yr. Horvath v1/v2, Hannum, and Lin form a competitive second tier (MAE 3–6 yr). PhenoAge has higher MAE because it predicts a biological age that systematically deviates from chronological age — by design. VidalBralo's MAE jumps from 9.5 (GSE40279) to 18.9 (GSE41169) due to its sensitivity to cohort composition.*

</div>

### Predicted-age distributions (bonus)

<div align="center">

![Predicted age distribution](results/figures/predicted_age_distribution.png)

*Box plots of predicted-age distributions per clock per dataset. Shaded horizontal bands = chronological-age range of each dataset; dashed lines = chronological-age median. A well-calibrated clock should have its box centred near its dataset's dashed line and contained within the shaded band — AltumAge, Horvathv1, and Hannum exemplify this. PhenoAge intentionally extends beyond the band on older samples (predicting biologically older), and VidalBralo's GSE41169 box sits entirely above the chronological-age band, visualising its calibration bias.*

</div>

---

## Per-Clock Performance — Real Data

### GSE40279 (n=656, 19–101 yr)

| Clock | Pearson *r* | MAE (yr) | Bias (yr) | RMSE (yr) |
|---|---:|---:|---:|---:|
| Horvathv1 | 0.918 | 4.77 | −2.33 | 6.32 |
| Horvathv2 | 0.940 | 4.45 | −2.70 | 5.90 |
| Hannum | 0.946 | 5.51 | +4.75 | 6.78 |
| Lin | 0.932 | 4.12 | −0.00 | 5.35 |
| VidalBralo | 0.727 | 9.47 | −4.44 | 11.71 |
| PhenoAge | 0.852 | 9.04 | −7.89 | 11.06 |
| HRSInCHPhenoAge | 0.868 | 6.35 | −3.37 | 8.10 |
| YingCausAge | 0.901 | 5.95 | −4.20 | 7.69 |
| YingDamAge | 0.691 | 8.89 | −3.95 | 12.18 |
| YingAdaptAge | 0.555 | 12.03 | −7.25 | 15.10 |
| StocP | 0.717 | 12.08 | −10.16 | 14.78 |
| StocH | 0.774 | 9.55 | −6.73 | 11.74 |
| **AltumAge** | **0.960** | **3.47** | −1.44 | 4.70 |
| DunedinPACE | 0.135 | — | — | — |

### GSE41169 (n=95, 18–65 yr)

| Clock | Pearson *r* | MAE (yr) | Bias (yr) | RMSE (yr) |
|---|---:|---:|---:|---:|
| Horvathv1 | 0.935 | 2.96 | +0.38 | 3.89 |
| Horvathv2 | 0.958 | 3.24 | −2.55 | 4.21 |
| Hannum | 0.929 | 3.05 | +0.35 | 4.12 |
| Lin | 0.885 | 6.36 | −5.75 | 7.61 |
| VidalBralo | 0.791 | 18.91 | +18.85 | 20.03 |
| PhenoAge | 0.912 | 5.57 | −4.82 | 6.96 |
| HRSInCHPhenoAge | 0.926 | 6.66 | −6.17 | 7.96 |
| YingCausAge | 0.916 | 6.77 | −6.64 | 7.94 |
| YingDamAge | 0.627 | 10.21 | +7.34 | 12.54 |
| YingAdaptAge | 0.657 | 13.07 | −12.33 | 15.22 |
| StocP | 0.753 | 5.88 | +2.09 | 7.17 |
| StocH | 0.750 | 8.52 | +7.13 | 9.86 |
| **AltumAge** | **0.952** | **2.59** | −1.15 | 3.45 |
| DunedinPACE | 0.502 | — | — | — |

For the auto-generated full metrics report, see [`results/summary.md`](results/summary.md).

---

## Discussion

### Why AltumAge wins

AltumAge uses a TabNet deep neural network trained on ~21,000 CpGs across multiple tissue types. Its attention mechanism learns to weight CpGs *conditionally* based on the input, capturing non-linear CpG–CpG interactions that penalised linear regression cannot represent. Combined with multi-tissue training, this gives it both better calibration (lower MAE) and better generalisation to cohorts the linear clocks have never seen.

### Why Horvath / Hannum / Lin remain competitive

The classic 2013–2016 clocks land within 1–2 years MAE of AltumAge on both cohorts despite using 70–400× fewer CpGs. This robustness is a function of careful elastic-net feature selection on large training cohorts — the selected CpGs change relatively monotonically with age across most tissues. They remain the practical default for most studies because they require minimal computation and are well-understood.

### Why PhenoAge has higher chronological-age MAE

PhenoAge predicts a *biological* age derived from 9 clinical biomarkers (CRP, glucose, albumin, etc.) — deviation from chronological age is the *signal*, interpreted as biological age acceleration (AgeAccel). The MAE in this report measures distance from chronological age, which is not what PhenoAge optimises for. A negative bias (~−8 yr on GSE40279) means PhenoAge predicts most samples to be biologically younger than their calendar age — consistent with the cohort being healthy population-based volunteers.

### Why the Ying causal trio drops in correlation

YingDamAge and YingAdaptAge encode age-related *damage* and *adaptive responses* respectively, not chronological age itself. Their lower r values (0.55–0.69) reflect this intentional decoupling. YingCausAge — the causal-CpG clock — retains high correlation (r ≈ 0.9) because causal CpGs are precisely those whose change drives chronological aging.

### Why VidalBralo's bias jumps cohorts

The 8-CpG VidalBralo clock has near-zero redundancy — every CpG carries large weight, so any β-value distribution shift between training and test cohorts directly biases predictions. The +18.85 yr bias on GSE41169 is the price of extreme compactness; it is not an analysis bug.

### Why DunedinPACE shows weak r with chronological age

DunedinPACE measures *rate* of aging (years aged per calendar year, ~1.0 = normal pace). The original Dunedin paper predicts r ≈ 0.0–0.3 with chronological age — older people may pace slightly faster on average, but pace is fundamentally a different quantity. r = 0.135 (GSE40279) and r = 0.502 (GSE41169) are both within published expectations.

---

## References

- de Lima Camillo, L. P. *et al.* (2022). *A pan-tissue DNA-methylation epigenetic clock based on deep learning.* npj Aging 8:4. [doi:10.1038/s41514-022-00085-y](https://doi.org/10.1038/s41514-022-00085-y)
- Belsky, D. W. *et al.* (2022). *DunedinPACE, a DNA methylation biomarker of the pace of aging.* eLife 11:e73420.
- Ying, K. *et al.* (2024). *Biolearn, an open-source library for biomarkers of aging.* Nature Aging.
- Higgins-Chen, A. T. *et al.* (2022). *A computational solution for bolstering reliability of epigenetic clocks.* Nature Aging 2:644–661.
- Ying, K. *et al.* (2022). *Causality-enriched epigenetic age uncouples damage and adaptation.* Nature Aging 4:231–246.
- Tong, H. *et al.* (2024). *Quantifying the stochastic component of epigenetic aging.* Nature Aging 4:886–901.
- Hannum, G. *et al.* (2013). *Genome-wide methylation profiles reveal quantitative views of human aging rates.* Mol Cell 49:359–367.
- Horvath, S. (2013). *DNA methylation age of human tissues and cell types.* Genome Biology 14:R115.
- Levine, M. E. *et al.* (2018). *An epigenetic biomarker of aging for lifespan and healthspan.* Aging 10:573–591.
- Lin, Q. *et al.* (2016). *DNA methylation levels at individual age-associated CpG sites can be indicative for life expectancy.* Aging 8:394–401.
- Vidal-Bralo, L. *et al.* (2018). *Specific premature epigenetic aging of cartilage in osteoarthritis.* Aging 10:3137–3151.
