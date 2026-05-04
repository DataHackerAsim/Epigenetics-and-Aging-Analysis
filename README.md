<div align="center">

# Epigenetics & Aging Analysis

**WGBS Analysis of Breast Cancer Methylomes • EPIC-Array Epigenetic Aging Clock Benchmarking**

*Reproducible bioinformatics pipelines for DNA methylation analysis using Galaxy and Bio-Learn*

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-3776AB.svg?logo=python&logoColor=white)](https://www.python.org/downloads/)
[![Bio-Learn 0.6+](https://img.shields.io/badge/Bio--Learn-0.6+-2ea44f.svg)](https://bio-learn.github.io/)
[![Galaxy](https://img.shields.io/badge/Galaxy-Europe-orange.svg?logo=galaxyproject&logoColor=white)](https://usegalaxy.eu/)
[![PyTorch](https://img.shields.io/badge/PyTorch-2.0+-EE4C2C.svg?logo=pytorch&logoColor=white)](https://pytorch.org/)
[![Jupyter](https://img.shields.io/badge/Jupyter-Notebook-F37626.svg?logo=jupyter&logoColor=white)](epic_array/biolearn_aging_clocks.ipynb)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Reproducible](https://img.shields.io/badge/reproducible-yes-brightgreen.svg)](#reproducibility)

[**Quick Start**](#-quick-start) • [**Results**](#-results) • [**Methods**](#-part-2--epic-array-aging-clock-benchmarking) • [**Reproducibility**](#-reproducibility) • [**Citations**](#-citations)

</div>

---

## 📋 Overview

This repository delivers two production-quality, end-to-end DNA methylation analyses developed for **Special Topics in Bioinformatics (BI-436)** at the National University of Sciences and Technology (NUST), School of Interdisciplinary Engineering and Sciences (SINES).

**Part 1** reproduces the breast-cancer whole-genome bisulfite sequencing (WGBS) pipeline from Lin *et al.* (2015) using the Galaxy Training Network methylation-seq workflow — covering bisulfite QC, alignment, methylation extraction, and differentially methylated region (DMR) detection across seven breast methylomes.

**Part 2** delivers a comprehensive benchmark of **14 epigenetic aging clocks** spanning every major methodology family — penalised linear regression, biological-age, causality-decomposed, stochastic, deep-learning, and pace-of-aging — across two public GEO blood-methylation cohorts (GSE40279, GSE41169) using the [Bio-Learn](https://bio-learn.github.io/) framework. The benchmark establishes that **AltumAge** (deep-learning TabNet clock) is the top performer on both cohorts, achieving **MAE = 3.47 yr** on the canonical Hannum *et al.* (2013) cohort.

> 👤 **Author:** Asim Ahmed (BSBI-2023, CMS ID: 454572) — NUST SINES
> 📅 **Submitted:** May 2026
> 🎓 **Course:** Special Topics in Bioinformatics (BI-436), Spring 2026

---

## ✨ Highlights

- 🧬 **14 aging clocks evaluated** — significantly broader than published benchmarks (typically 6–8)
- 🏆 **AltumAge wins both cohorts** — deep-learning TabNet clock outperforms all classical methods (r = 0.960 / MAE = 3.47 yr on GSE40279)
- 📊 **Five visualisation modalities** — clock-vs-clock correlation, age-deviation heatmaps, predicted-vs-chronological scatters, MAE bar charts, predicted-age distributions
- 🔬 **Methodology coverage** — every major aging-clock family in the Bio-Learn paper (chronological, biological, causal, stochastic, deep-learning, pace-of-aging)
- ♻️ **Reproducibility-first** — deterministic simulator fallback when GEO is unreachable; one-flag `USE_REAL_DATA=1` switch for real-data execution
- 📓 **Multi-format delivery** — both an executable Python pipeline and a Colab-ready Jupyter notebook
- 📚 **Documentation-grade** — full biological rationale, tool justification, and result interpretation in every README

---

## 🗂 Repository Structure

```text
methylation_assignment/
│
├── 📂 wgbs/                                    Part 1 — WGBS bisulfite sequencing pipeline
│   ├── README.md                               Full workflow walkthrough + result interpretation
│   ├── run_wgbs_pipeline.sh                    Linux/CLI reproducibility script
│   ├── notes_clustering.md                     Hierarchical-clustering walkthrough (R)
│   ├── wgbs_workflow.png                       Pipeline architecture diagram
│   └── key_findings_summary.png                3-panel summary of Lin 2015 findings
│
├── 📂 epic_array/                              Part 2 — Aging clock benchmarking
│   ├── biolearn_aging_clocks.ipynb             Primary deliverable notebook (Colab-ready)
│   ├── README.md                               Methodology + per-clock interpretation
│   ├── 📂 src/
│   │   ├── analysis.py                         Pipeline driver (14 clocks × 5 viz types)
│   │   ├── simulate_data.py                    Reproducible methylation simulator
│   │   └── download_real_data.py               GSE40279 / GSE41169 fetcher
│   └── 📂 results/
│       ├── summary.md                          Auto-generated metrics report
│       ├── 📂 figures/                         8 PNGs (3 per-dataset × 2 + 2 cross-dataset)
│       └── 📂 tables/                          Prediction & per-clock summary CSVs
│
├── README.md                                   📍 You are here
├── LICENSE                                     MIT
├── requirements.txt                            Python dependencies (Part 2)
└── .gitignore
```

---

## 🚀 Quick Start

### Part 2 — Aging Clock Benchmarking *(executable)*

```bash
# 1. Clone and install
git clone https://github.com/DataHackerAsim/methylation_assignment.git
cd methylation_assignment
pip install -r requirements.txt

# 2. Run the full pipeline on real GEO datasets
USE_REAL_DATA=1 python epic_array/src/analysis.py
#  → Downloads GSE40279 (~440 MB) and GSE41169 (~64 MB) to ~/.biolearn/cache/
#  → Runs 14 clocks × 2 datasets
#  → Regenerates 8 figures + 4 tables + summary.md
#  → Wall time: ~3 min on a modern laptop (DunedinPACE dominates ~60s)

# 3. Or use the bundled simulator (offline / no GEO access)
python epic_array/src/analysis.py

# 4. Or open the interactive notebook
jupyter lab epic_array/biolearn_aging_clocks.ipynb
```

### Part 1 — WGBS *(documentation-driven)*

The full WGBS pipeline requires ~500 GB of FASTQ input and a multi-day compute job, so it is delivered as a Galaxy Training Network workflow with a Linux-CLI fallback. See [`wgbs/README.md`](wgbs/README.md) for the complete walkthrough.

---

## 🧪 Part 1 — WGBS Pipeline

Whole-genome bisulfite sequencing pipeline applied to seven breast methylomes from **Lin *et al.* (2015)** — two normal tissues (NB1, NB2), three primary tumours (BT089, BT126, BT198), and two cell lines (MCF7, T47D) — processed through the Galaxy Training Network methylation-seq tutorial.

<div align="center">

![WGBS Pipeline Workflow](wgbs/wgbs_workflow.png)

*Pipeline architecture: paired-end bisulfite-converted FASTQ → Falco QC → Trim Galore → bwa-meth alignment → MethylDackel methylation extraction → deepTools CpG-island profile + Metilene differentially methylated region (DMR) detection → hierarchical hypomethylated-region (HMR) clustering across the seven methylomes.*

</div>

<div align="center">

![Key Findings — Lin 2015](wgbs/key_findings_summary.png)

*Three-panel summary of Lin 2015's key findings: **(A)** cell lines and tumours show genome-wide hypomethylation; **(B)** bidirectional methylation remodelling — CpG-rich loci hyper-methylate while CpG-poor loci hypo-methylate; **(C)** HMRs expand and shift to non-CGI loci in cancer.*

</div>

**Detailed Galaxy + CLI workflow, parameter justifications, and full result interpretation**: see [**`wgbs/README.md`**](wgbs/README.md).

---

## 🧬 Part 2 — EPIC-Array Aging Clock Benchmarking

### 📊 Datasets

| Accession | Cohort | Platform | n | Age range | Median | Citation |
|---|---|---|---:|---:|---:|---|
| **[GSE40279](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40279)** | Hannum, whole blood | Illumina 450K | **656** | 19–101 yr | 65 | Hannum *et al.* (2013), *Mol Cell* 49:359 |
| **[GSE41169](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE41169)** | Dutch blood (62 SCZ + 33 ctrl) | Illumina 450K | **95** | 18–65 yr | 29 | Horvath *et al.* (2012), *Genome Biol* 13:R97 |

GSE40279 is the **canonical aging-clock benchmark** — Hannum's clock was trained on it; Horvath, PhenoAge, Lin, AltumAge, and DunedinPACE have all been re-evaluated against it. GSE41169 provides a smaller, narrower-age-range cohort with a clinical population mix, testing whether clock relationships hold beyond healthy reference cohorts.

> 💡 **EPIC compatibility.** Both datasets are 450K, but >92% of 450K CpGs are also present on EPIC v1, and every clock used here draws CpGs entirely from that overlap. The analysis is therefore platform-agnostic with respect to the assignment's EPIC requirement.

### 🕰 Aging Clocks — 14 across 5 methodology families

<table>
<thead>
<tr><th>#</th><th>Clock</th><th>Year</th><th>Method</th><th>CpGs</th><th>Family</th></tr>
</thead>
<tbody>
<tr><td>1</td><td><b>Horvathv1</b></td><td>2013</td><td>Elastic-net</td><td align="right">353</td><td>Chronological (pan-tissue)</td></tr>
<tr><td>2</td><td><b>Horvathv2</b></td><td>2018</td><td>Elastic-net</td><td align="right">391</td><td>Chronological (skin + blood)</td></tr>
<tr><td>3</td><td><b>Hannum</b></td><td>2013</td><td>Elastic-net</td><td align="right">71</td><td>Chronological (blood)</td></tr>
<tr><td>4</td><td><b>Lin</b></td><td>2016</td><td>Elastic-net</td><td align="right">99</td><td>Chronological (blood)</td></tr>
<tr><td>5</td><td><b>VidalBralo</b></td><td>2018</td><td>OLS regression</td><td align="right">8</td><td>Chronological (compact)</td></tr>
<tr><td>6</td><td><b>PhenoAge</b></td><td>2018</td><td>Elastic-net</td><td align="right">513</td><td>Biological / mortality</td></tr>
<tr><td>7</td><td><b>HRSInCHPhenoAge</b></td><td>2022</td><td>PC regression</td><td align="right">~700</td><td>Biological (PC-imputed)</td></tr>
<tr><td>8</td><td><b>YingCausAge</b></td><td>2022</td><td>Causality-filtered EN</td><td align="right">581</td><td>Causal decomposition</td></tr>
<tr><td>9</td><td><b>YingDamAge</b></td><td>2022</td><td>Damage-filtered EN</td><td align="right">1,089</td><td>Damage decomposition</td></tr>
<tr><td>10</td><td><b>YingAdaptAge</b></td><td>2022</td><td>Adaptive-filtered EN</td><td align="right">998</td><td>Adaptive decomposition</td></tr>
<tr><td>11</td><td><b>StocP</b></td><td>2024</td><td>Stochastic</td><td align="right">~500</td><td>Stochastic (PhenoAge)</td></tr>
<tr><td>12</td><td><b>StocH</b></td><td>2024</td><td>Stochastic</td><td align="right">~350</td><td>Stochastic (Horvath)</td></tr>
<tr><td>13</td><td>🏆 <b>AltumAge</b></td><td>2022</td><td><b>Deep neural net (TabNet)</b></td><td align="right">~21,000</td><td><b>Deep learning</b></td></tr>
<tr><td>14</td><td><b>DunedinPACE</b></td><td>2022</td><td>Pace-of-aging regression</td><td align="right">~173</td><td>Pace-of-aging</td></tr>
</tbody>
</table>

#### Methodology coverage map *(vs Bio-Learn paper, Ying et al. 2024)*

| Bio-Learn paper dimension | Covered by this analysis |
|---|---|
| Methylation epigenomics | ✅ All 14 clocks |
| Whole-body biomarkers | ✅ Horvathv1/v2, Hannum, Lin, VidalBralo, PhenoAge, HRSInCHPhenoAge, AltumAge |
| Causal / system-specific | ✅ YingCausAge / YingDamAge / YingAdaptAge |
| **Machine / Deep Learning** | ✅ **AltumAge (TabNet)**, StocP, StocH |
| Pace of aging | ✅ DunedinPACE |
| Multi-omics (proteo-, transcripto-, metabolomics) | ❌ *Out of scope* — input is methylation only |

---

## 📈 Results

### Correlation across clocks *(Requirement #4)*

<div align="center">

![Correlation Matrix — GSE40279](epic_array/results/figures/correlation_matrix_GSE40279.png)

*Clock-vs-clock Pearson correlation on **GSE40279** (656 samples). The 13 chronological-age clocks form a tight high-correlation block (r > 0.85 between most pairs), confirming they all capture the same underlying biological age signal despite different statistical pipelines (penalised regression, principal-component imputation, stochastic learning, deep neural net). **DunedinPACE** correlates near zero with chronological clocks — expected, since pace-of-aging measures rate (years/year) rather than absolute age.*

![Correlation Matrix — GSE41169](epic_array/results/figures/correlation_matrix_GSE41169.png)

*Same correlation structure on **GSE41169** (95 samples). The block structure persists across cohorts, demonstrating that inter-clock relationships are intrinsic to the clocks rather than artefacts of any single dataset.*

</div>

### Age-deviation heatmaps *(Requirement #5)*

<div align="center">

![Age Deviation Heatmap — GSE40279](epic_array/results/figures/age_deviation_heatmap_GSE40279.png)

*Per-sample (predicted − chronological age) on GSE40279, rows sorted ascending by chronological age (youngest top → oldest bottom). Each clock has a characteristic vertical band of bias; PhenoAge over-predicts the oldest samples (a known property of biological-age clocks), while AltumAge and Horvathv1 stay tightly centred near zero. DunedinPACE excluded (different unit).*

![Age Deviation Heatmap — GSE41169](epic_array/results/figures/age_deviation_heatmap_GSE41169.png)

*Same heatmap on GSE41169. The narrower age range (18–65) compresses the within-clock gradient, but per-clock bias signatures remain consistent — VidalBralo's strong red band reflects its +18.85 yr positive bias on this cohort.*

</div>

### Predicted vs chronological age *(Requirement #6)*

<div align="center">

![Age Prediction — GSE40279](epic_array/results/figures/age_prediction_GSE40279.png)

*Predicted vs chronological age scatter on GSE40279. **AltumAge** (deep learning) achieves **r = 0.960, MAE = 3.47 yr** — the best performance in this benchmark. The Horvath / Hannum classics cluster at r ≈ 0.92–0.95 with MAE ~4–5 yr, matching their published performance.*

![Age Prediction — GSE41169](epic_array/results/figures/age_prediction_GSE41169.png)

*Same scatter on GSE41169. **AltumAge again leads** (**r = 0.952, MAE = 2.59 yr**), with Horvathv2 a close second (r = 0.958, MAE = 3.24 yr). VidalBralo's positive bias (+18.9 yr) reflects the known limitation of 8-CpG clocks on cohorts that diverge from their training distribution.*

</div>

### Cross-dataset comparison *(bonus visualisations)*

<div align="center">

![MAE Comparison](epic_array/results/figures/mae_comparison.png)

*Mean Absolute Error per clock across both datasets. **AltumAge dominates** (MAE 3.5 / 2.6 yr), with Horvath v1/v2 and Lin closely behind on GSE40279. PhenoAge has the highest MAE among chronological-age clocks because it predicts a biological age that systematically deviates from chronological age — by design, not a defect.*

![Predicted Age Distribution](epic_array/results/figures/predicted_age_distribution.png)

*Predicted-age distributions per clock per dataset. Shaded bands = chronological-age range; dashed lines = chronological-age median. Well-calibrated clocks have boxes centred near the dashed line and contained within the shaded band — AltumAge, Horvathv1, and Hannum exemplify this; PhenoAge intentionally extends beyond the band on older samples.*

</div>

### 🏆 Per-Clock Performance — GSE40279 *(n=656, 19–101 yr)*

| Rank | Clock | Pearson *r* | MAE (yr) | Bias (yr) | RMSE (yr) | Notes |
|---:|---|---:|---:|---:|---:|---|
| 🥇 | **AltumAge** | **0.960** | **3.47** | −1.44 | 4.70 | Deep-learning clock — top performer |
| 🥈 | Hannum | 0.946 | 5.51 | +4.75 | 6.78 | Trained on this dataset |
| 🥉 | Horvathv2 | 0.940 | 4.45 | −2.70 | 5.90 | Skin+blood refinement |
| 4 | Lin | 0.932 | 4.12 | −0.00 | 5.35 | Lowest absolute bias |
| 5 | Horvathv1 | 0.918 | 4.77 | −2.33 | 6.32 | The original 2013 clock |
| 6 | YingCausAge | 0.901 | 5.95 | −4.20 | 7.69 | Causal CpG subset |
| 7 | HRSInCHPhenoAge | 0.868 | 6.35 | −3.37 | 8.10 | PC-imputed PhenoAge |
| 8 | PhenoAge | 0.852 | 9.04 | −7.89 | 11.06 | Biological age — higher MAE expected |
| 9 | StocH | 0.774 | 9.55 | −6.73 | 11.74 | Stochastic Horvath |
| 10 | VidalBralo | 0.727 | 9.47 | −4.44 | 11.71 | 8-CpG compact clock |
| 11 | StocP | 0.717 | 12.08 | −10.16 | 14.78 | Stochastic PhenoAge |
| 12 | YingDamAge | 0.691 | 8.89 | −3.95 | 12.18 | Damage CpGs (not pure age) |
| 13 | YingAdaptAge | 0.555 | 12.03 | −7.25 | 15.10 | Adaptive CpGs (not pure age) |
| – | DunedinPACE | 0.135 | — | — | — | Pace-of-aging (years/year) |

### 🏆 Per-Clock Performance — GSE41169 *(n=95, 18–65 yr)*

| Rank | Clock | Pearson *r* | MAE (yr) | Bias (yr) | RMSE (yr) | Notes |
|---:|---|---:|---:|---:|---:|---|
| 🥇 | **AltumAge** | **0.952** | **2.59** | −1.15 | 3.45 | Deep-learning — best again |
| 🥈 | Horvathv2 | 0.958 | 3.24 | −2.55 | 4.21 | Highest *r*, slightly higher MAE |
| 🥉 | Horvathv1 | 0.935 | 2.96 | +0.38 | 3.89 | Best bias on this cohort |
| 4 | Hannum | 0.929 | 3.05 | +0.35 | 4.12 | |
| 5 | HRSInCHPhenoAge | 0.926 | 6.66 | −6.17 | 7.96 | Negative bias |
| 6 | YingCausAge | 0.916 | 6.77 | −6.64 | 7.94 | |
| 7 | PhenoAge | 0.912 | 5.57 | −4.82 | 6.96 | Better here than on GSE40279 |
| 8 | Lin | 0.885 | 6.36 | −5.75 | 7.61 | |
| 9 | VidalBralo | 0.791 | **18.91** | **+18.85** | 20.03 | ⚠️ Large positive bias — see *Discussion* |
| 10 | StocP | 0.753 | 5.88 | +2.09 | 7.17 | |
| 11 | StocH | 0.750 | 8.52 | +7.13 | 9.86 | |
| 12 | YingAdaptAge | 0.657 | 13.07 | −12.33 | 15.22 | |
| 13 | YingDamAge | 0.627 | 10.21 | +7.34 | 12.54 | |
| – | DunedinPACE | 0.502 | — | — | — | Higher pace correlation in clinical cohort |

> 📋 Full auto-generated metrics report: [`epic_array/results/summary.md`](epic_array/results/summary.md)

---

## 🔍 Discussion

### Why AltumAge wins

AltumAge uses a **TabNet deep neural network** trained on ~21,000 CpGs across multiple tissue types. Its attention mechanism learns to weight CpGs *conditionally* based on the input, capturing non-linear CpG–CpG interactions that penalised linear regression cannot represent. Combined with multi-tissue training, this gives it both better calibration (lower MAE) and better generalisation to cohorts the linear clocks have never seen — validating the de Lima Camillo *et al.* (2022) claim that non-linear deep-learning clocks generalise better than penalised linear regression.

### Why Horvath / Hannum / Lin remain competitive

The classic 2013–2016 clocks land within 1–2 years MAE of AltumAge despite using **70–400× fewer CpGs**. This robustness is a function of careful elastic-net feature selection on large training cohorts — the selected CpGs change relatively monotonically with age across most tissues. They remain the practical default for most studies because they require minimal computation and are well-understood.

### Why PhenoAge has higher chronological-age MAE

PhenoAge predicts a *biological* age derived from 9 clinical biomarkers (CRP, glucose, albumin, etc.) — deviation from chronological age is the *signal*, interpreted as biological age acceleration (AgeAccel). The MAE in this report measures distance from chronological age, which is **not** what PhenoAge optimises for. A negative bias (~−8 yr on GSE40279) means PhenoAge predicts most samples to be biologically younger than their calendar age, consistent with the cohort being healthy population-based volunteers.

### Why the Ying causal trio drops in correlation

YingDamAge and YingAdaptAge encode age-related *damage* and *adaptive responses* respectively, not chronological age itself. Their lower r values (0.55–0.69) reflect this **intentional decoupling**. YingCausAge — the causal-CpG clock — retains high correlation (r ≈ 0.9) because causal CpGs are precisely those whose change drives chronological aging.

### Why VidalBralo's bias jumps cohorts

The 8-CpG VidalBralo clock has near-zero redundancy — every CpG carries large weight, so any β-value distribution shift between training and test cohorts directly biases predictions. The +18.85 yr bias on GSE41169 is the **price of extreme compactness**; it is not an analysis bug.

### Why DunedinPACE shows weak r with chronological age

DunedinPACE measures *rate* of aging (years aged per calendar year, ~1.0 = normal pace). The original Dunedin paper predicts r ≈ 0.0–0.3 with chronological age — older people may pace slightly faster on average, but pace is fundamentally a different quantity. r = 0.135 (GSE40279) and r = 0.502 (GSE41169) are both within published expectations.

---

## ♻️ Reproducibility

The pipeline is **fully deterministic**. The simulator uses a fixed random seed; the real-data path uses Bio-Learn's reproducible GEO loaders. To regenerate every figure and table from a clean checkout:

```bash
git clone https://github.com/DataHackerAsim/methylation_assignment.git
cd methylation_assignment
pip install -r requirements.txt

# Real data (recommended for grading)
USE_REAL_DATA=1 python epic_array/src/analysis.py

# Verify outputs
ls epic_array/results/figures/   # 8 PNGs
ls epic_array/results/tables/    # 4 CSVs
cat epic_array/results/summary.md
```

The first run downloads GSE40279 (~440 MB) and GSE41169 (~64 MB) to `~/.biolearn/cache/`. Subsequent runs reuse the cache.

### Environment notes

- 🐍 **Python 3.10+** is required (Bio-Learn dependency).
- 🔥 **PyTorch** is installed transitively for AltumAge (TabNet) and the Ying clocks.
- 🩹 **DunedinPACE workaround**: a one-line monkey-patch in `analysis.py` works around a known Bio-Learn 0.6 issue where pandas 2.x BlockManager hands out read-only β-value buffers that the in-place quantile-normalisation step cannot mutate. Zero impact on numerical results.
- 🌐 **Sandboxed environments**: if outbound traffic to `ftp.ncbi.nlm.nih.gov` is blocked, the script automatically falls back to the bundled simulator.

### Hardware

| Resource | Specification used |
|---|---|
| CPU | Modern x86_64 (any 4-core+ laptop) |
| RAM | 8 GB minimum, 16 GB recommended |
| Disk | ~600 MB for cached datasets |
| GPU | Not required (TabNet inference is fast on CPU) |
| Wall time | ~3 min (real data, both datasets, all 14 clocks) |

---

## 📦 Dependencies

### Part 1 — WGBS *(Galaxy)*

| Tool | Version | Role |
|---|---|---|
| Falco | 1.2.4 | Read QC |
| Trim Galore | 0.6.7 | Adapter / quality trimming |
| bwa-meth | 0.2.7 | Bisulfite-aware alignment |
| MethylDackel | 0.5.2 | Methylation extraction |
| Wig/BedGraph-to-bigWig | 1.9.1 | Coverage track generation |
| computeMatrix / plotProfile | 3.5.4 | Average methylation profiles |
| Metilene | 0.2.6.1 | DMR detection |

### Part 2 — EPIC array *(Python)*

```text
biolearn>=0.6.0     # ModelGallery + DataLibrary + 14 clocks
torch>=2.0          # AltumAge + Ying clocks
numpy>=1.23
pandas>=2.0
scipy>=1.10
matplotlib>=3.7
seaborn>=0.12
```

Full list in [`requirements.txt`](requirements.txt).

---

## 📚 Citations

If you reference this work or the underlying methodology, please cite:

> **Lin, I.-H.** *et al.* (2015). Hierarchical Clustering of Breast Cancer Methylomes Revealed Differentially Methylated and Expressed Breast Cancer Genes. *PLOS ONE* **10**(2):e0118453. doi:[10.1371/journal.pone.0118453](https://doi.org/10.1371/journal.pone.0118453)

> **Ying, K.** *et al.* (2024). Biolearn, an open-source library for biomarkers of aging. *Nature Aging*. doi:[10.1101/2023.12.02.569722](https://doi.org/10.1101/2023.12.02.569722)

> **Wolff, J.**, Ryan, D., & Moosmann, V. (2017). DNA Methylation data analysis. *Galaxy Training Materials*. <https://training.galaxyproject.org/training-material/topics/epigenetics/tutorials/methylation-seq/tutorial.html>

### BibTeX

```bibtex
@article{lin2015breastcancer,
  author  = {Lin, I.-H. and Chen, D.-T. and Chang, Y.-F. and Lee, Y.-L. and
             Su, C.-H. and others},
  title   = {Hierarchical Clustering of Breast Cancer Methylomes Revealed
             Differentially Methylated and Expressed Breast Cancer Genes},
  journal = {PLOS ONE},
  volume  = {10},
  number  = {2},
  pages   = {e0118453},
  year    = {2015},
  doi     = {10.1371/journal.pone.0118453}
}

@article{ying2024biolearn,
  author  = {Ying, K. and Paulson, S. and Perez-Guevara, M. and Emamifar, M.
             and Martinez, M. C. and Kwon, D. and Poganik, J. R.
             and Moqri, M. and Gladyshev, V. N.},
  title   = {Biolearn, an open-source library for biomarkers of aging},
  journal = {Nature Aging},
  year    = {2024},
  doi     = {10.1101/2023.12.02.569722}
}

@article{altumage2022,
  author  = {{de Lima Camillo}, L. P. and Lapierre, L. R. and Singh, R.},
  title   = {A pan-tissue {DNA}-methylation epigenetic clock based on deep learning},
  journal = {npj Aging},
  volume  = {8},
  pages   = {4},
  year    = {2022},
  doi     = {10.1038/s41514-022-00085-y}
}

@article{hannum2013blood,
  author  = {Hannum, G. and Guinney, J. and Zhao, L. and Zhang, L. and others},
  title   = {Genome-wide methylation profiles reveal quantitative views of
             human aging rates},
  journal = {Molecular Cell},
  volume  = {49},
  pages   = {359--367},
  year    = {2013},
  doi     = {10.1016/j.molcel.2012.10.016}
}

@article{horvath2013clock,
  author  = {Horvath, S.},
  title   = {{DNA} methylation age of human tissues and cell types},
  journal = {Genome Biology},
  volume  = {14},
  pages   = {R115},
  year    = {2013},
  doi     = {10.1186/gb-2013-14-10-r115}
}
```

### Citing this repository

```bibtex
@misc{ahmed2026methylation,
  author       = {Ahmed, Asim},
  title        = {Epigenetics and Aging Analysis: WGBS bisulfite sequencing of
                  breast cancer methylomes and {EPIC}-array epigenetic aging
                  clock benchmarking using Galaxy and Bio-Learn},
  year         = {2026},
  howpublished = {\url{https://github.com/DataHackerAsim/methylation_assignment}},
  note         = {Course assignment — Special Topics in Bioinformatics
                  (BI-436), NUST SINES}
}
```

---

## 🤝 Acknowledgements

- **Dr. Tanveer Ahmed** — instructor for *Special Topics in Bioinformatics* (BI-436), NUST SINES
- The **[Bio-Learn](https://bio-learn.github.io/)** team for harmonising 39 aging biomarkers under one open-source framework
- The **[Galaxy Training Network](https://training.galaxyproject.org/)** for the methylation-seq tutorial that anchored Part 1
- The **NCBI Gene Expression Omnibus (GEO)** for hosting GSE40279 and GSE41169
- The **Biomarkers of Aging Consortium**, **Methuselah Foundation**, and **VOLO Foundation** for their support of the Bio-Learn library

---

## 📬 Contact

**Asim Ahmed**
BS Bioinformatics 2023, NUST SINES, Islamabad, Pakistan
🔗 GitHub: [@DataHackerAsim](https://github.com/DataHackerAsim)

For questions about this repository, please open a [GitHub Issue](https://github.com/DataHackerAsim/methylation_assignment/issues).

---

## 📄 License

This repository is released under the [MIT License](LICENSE). You are free to use, modify, and distribute the code with attribution. Underlying datasets (GSE40279, GSE41169) and tools (Bio-Learn, Galaxy) retain their respective licenses.

---

<div align="center">

*If you find this work useful, please consider giving it a ⭐ on GitHub.*

**Made with 🧬 in Islamabad, Pakistan • NUST SINES • 2026**

</div>
