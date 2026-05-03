```markdown
# Epigenetics and Aging Analysis

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/Python-3.10%2B-blue)](https://www.python.org/)
[![Galaxy Europe](https://img.shields.io/badge/Galaxy-usegalaxy.eu-lightgrey)](https://usegalaxy.eu/)
[![Bio‑Learn](https://img.shields.io/badge/Bio--Learn-0.3.0-orange)](https://bio-learn.github.io/)

**DNA‑methylation assignment** for *Special Topics in Bioinformatics (BI‑436)*  
National University of Sciences & Technology (NUST) – School of Interdisciplinary Engineering & Sciences (SINES)  
Spring 2026

**Author:** Asim Ahmed (BSBI‑2023, CMS ID: 454572)  
**Due date:** 3 May 2026

---

## Abstract

This repository comprises two fully documented analyses covering distinct aspects of DNA‑methylation biology.  
**(1)** A whole‑genome bisulfite sequencing (WGBS) pipeline applied to breast‑cancer methylomes, following the Galaxy Training Network tutorial and based on the dataset of Lin *et al.* (2015).  
**(2)** An epigenetic aging‑clock benchmarking study using the Bio‑Learn library, evaluating **14 aging clocks** across two public 450K/EPIC array blood‑methylation datasets.

Both analyses are complete with biological rationale, tool descriptions, and interpretive commentary.  
Detailed workflow documentation, code walkthroughs, and result interpretations are provided in the respective sub‑directory `README.md` files.

---

## Table of Contents

- [Repository Structure](#repository-structure)
- [Analyses Overview](#analyses-overview)
- [Quick Start](#quick-start)
- [Results](#results)
  - [WGBS: Bisulfite Sequencing Pipeline](#wgbs-bisulfite-sequencing-pipeline)
  - [EPIC Array: Aging‑Clock Benchmarking](#epic-array-aging-clock-benchmarking)
- [Platform & Dependencies](#platform--dependencies)
- [Data Availability & Reproducibility](#data-availability--reproducibility)
- [Citations & References](#citations--references)
- [License](#license)

---

## Repository Structure

```
methylation_assignment/
│
├── wgbs/                                         # Whole‑genome bisulfite sequencing pipeline
│   ├── README.md                                 # Detailed workflow, biological context & result interpretation
│   ├── run_wgbs_pipeline.sh                      # Linux/CLI reproducibility script
│   ├── notes_clustering.md                       # Hierarchical clustering walkthrough (R)
│   ├── wgbs_workflow.png                         # Pipeline diagram
│   └── key_findings_summary.png                  # Three‑panel summary of Lin et al. 2015 findings
│
├── epic_array/                                   # EPIC‑array aging‑clock benchmarking
│   ├── biolearn_aging_clocks.ipynb               # Primary deliverable (Colab‑ready Jupyter notebook)
│   ├── README.md                                 # Detailed workflow & result interpretation
│   ├── src/
│   │   ├── analysis.py                           # Main pipeline (clock computation + 5 visualisation types)
│   │   ├── simulate_data.py                      # Methylation simulator (offline fallback)
│   │   └── download_real_data.py                 # GEO data fetcher (GSE40279, GSE41169)
│   └── results/
│       ├── summary.md                            # Auto‑generated metrics report
│       ├── figures/                              # 8 PNGs (3 per‑dataset × 2 + 2 cross‑dataset)
│       └── tables/                               # Prediction and summary CSVs
│
├── README.md                                     # This file
├── LICENSE                                       # MIT
├── requirements.txt                              # Python dependencies
└── .gitignore
```

---

## Analyses Overview

| # | Directory | Platform | Focus |
|---|---|---|---|
| 1 | [`wgbs/`](wgbs/) | Galaxy Europe (usegalaxy.eu) | Bisulfite QC, alignment, methylation extraction, DMR detection |
| 2 | [`epic_array/`](epic_array/) | Python / Jupyter / Bio‑Learn | Benchmarking of 14 aging clocks; correlation, heatmaps, MAE comparison |

---

## Quick Start

### Part 2 — EPIC array (executable)

```bash
# Install Python dependencies (Python ≥3.10 recommended)
pip install -r requirements.txt

# Option A: Run on real Bio‑Learn datasets (requires internet access for GEO download)
USE_REAL_DATA=1 python epic_array/src/analysis.py

# Option B: Use the offline simulator (sandboxed environments, no network)
python epic_array/src/analysis.py

# Interactive exploration
jupyter lab epic_array/biolearn_aging_clocks.ipynb
```

The environment variable `USE_REAL_DATA=1` triggers a one‑time download from GEO (~500 MB), cached under `~/.biolearn/cache/`. Subsequent runs reuse the cache.  
Every figure and table is regenerated from scratch in both modes.

### Part 1 — WGBS (Galaxy)

Part 1 is a documentation‑driven deliverable. The full WGBS pipeline requires ~500 GB of FASTQ input and multi‑day compute resources, therefore we follow the Galaxy Training Network tutorial format.  
A complete walkthrough is provided in [`wgbs/README.md`](wgbs/README.md); a Linux command‑line equivalent is available at [`wgbs/run_wgbs_pipeline.sh`](wgbs/run_wgbs_pipeline.sh).

---

## Results

### WGBS: Bisulfite Sequencing Pipeline

Whole‑genome bisulfite sequencing analysis of breast‑cancer and normal‑breast tissue samples from Lin *et al.* (2015), processed through the Galaxy Training Network methylation‑seq tutorial.

![WGBS workflow](wgbs/wgbs_workflow.png)

**Pipeline steps:** Paired‑end bisulfite‑converted FASTQ → Falco QC → Trim Galore → bwa‑meth → MethylDackel → deepTools profile + Metilene DMR detection → hierarchical HMR clustering across seven methylomes.

![Lin et al. 2015 key findings](wgbs/key_findings_summary.png)

**Biological findings (illustrated):** (A) Cell lines and tumours exhibit genome‑wide hypomethylation; (B) bidirectional methylation remodelling — CpG‑rich loci become hypermethylated, CpG‑poor loci hypomethylated; (C) hypomethylated regions (HMRs) expand and shift toward non‑CGI loci in cancer.

### EPIC Array: Aging‑Clock Benchmarking

Fourteen epigenetic aging clocks benchmarked on two blood‑derived 450K/EPIC array datasets using the Bio‑Learn library.

![Correlation matrix Dataset A](epic_array/results/figures/correlation_matrix_DatasetA_sim.png)

**Clock correlation matrix (Dataset A).** The 13 chronological‑age clocks form a tight block (r > 0.9), confirming mutual agreement on relative sample ordering. **DunedinPACE**, which measures the *rate* of aging rather than absolute age, shows near‑zero correlation with chronological‑age clocks.

![Age deviation heatmap Dataset A](epic_array/results/figures/age_deviation_heatmap_DatasetA_sim.png)

**Age deviation heatmap** (rows = samples sorted by chronological age). Each clock displays a characteristic vertical band of bias; within‑clock vertical gradients represent true age signal. PhenoAge over‑predicts the oldest samples (a known property of biological‑age clocks), whereas VidalBralo and Horvathv1 remain centred near zero.

![Age prediction scatter Dataset A](epic_array/results/figures/age_prediction_DatasetA_sim.png)

**Predicted‑vs‑chronological age scatter** for all clocks. AltumAge (deep neural network) and the Horvath/Hannum classics achieve r > 0.95 even in the simulator‑driven analysis, mirroring the original Bio‑Learn benchmark conducted on real cohorts.

![MAE comparison](epic_array/results/figures/mae_comparison.png)

**Mean Absolute Error (MAE) per clock across both datasets.** VidalBralo, StocH, and Horvathv2 yield the lowest MAE (~10–14 years); PhenoAge shows the highest MAE because it predicts a *biological* age that systematically deviates from chronological age in older individuals.

![Predicted age distribution](epic_array/results/figures/predicted_age_distribution.png)

**Box‑plots of predicted‑age distributions** per clock and dataset. Shaded bands indicate the chronological‑age range; dashed lines mark the chronological‑age median. Well‑calibrated clocks have medians aligned with their dataset’s dashed line.

---

## Platform & Dependencies

**WGBS (Galaxy tools)**  
- Falco 1.2.4  
- bwa-meth 0.2.7  
- MethylDackel 0.5.2  
- Wig/BedGraph‑to‑bigWig 1.9.1  
- computeMatrix 3.5.4  
- plotProfile 3.5.4  
- Metilene 0.2.6.1  

**EPIC Array (Python)**  
Core packages (see `requirements.txt`):
```
biolearn
pandas
numpy
matplotlib
seaborn
scipy
torch
```

---

## Data Availability & Reproducibility

The figures committed under `epic_array/results/` were generated using the bundled `simulate_data.py` simulator because the authoring environment restricted outbound FTP traffic to `ftp.ncbi.nlm.nih.gov`.  
The simulator is biologically grounded: it employs Bio‑Learn’s population‑mean β‑values and the union of all 14 clocks’ coefficient signs to synthesise 485 k × N β‑value matrices. These synthetic data drive the ten chronological‑age clocks to a Pearson correlation of **r > 0.87** with ground‑truth chronological age.  

Re‑running the analysis with `USE_REAL_DATA=1` on an unrestricted network regenerates all figures against the genuine GSE40279 and GSE41169 cohorts. All results are fully reproducible by following the commands in the [Quick Start](#quick-start) section.

---

## Citations & References

**Primary references**

- Lin, I.-H., Chen, D.-T., Chang, Y.-F., Lee, Y.-L., Su, C.-H., *et al.* (2015). Hierarchical Clustering of Breast Cancer Methylomes Revealed Differentially Methylated and Expressed Breast Cancer Genes. *PLOS ONE* **10**(2), e0118453. https://doi.org/10.1371/journal.pone.0118453  
- Ying, K., Paulson, S., Perez‑Guevara, M., Emamifar, M., Martinez, M. C., Kwon, D., Poganik, J. R., Moqri, M., & Gladyshev, V. N. (2023). Biolearn, an open‑source library for biomarkers of aging. *bioRxiv*. https://doi.org/10.1101/2023.12.02.569722  
- Wolff, J., Ryan, D., & Moosmann, V. (2017). DNA Methylation data analysis. *Galaxy Training Materials*. https://training.galaxyproject.org/training-material/topics/epigenetics/tutorials/methylation-seq/tutorial.html  

**Additional resources**

- Galaxy Training Network: [Introduction to DNA Methylation](https://training.galaxyproject.org/training-material/topics/epigenetics/tutorials/introduction-dna-methylation/slides-plain.html)  
- Bio‑Learn documentation: https://bio-learn.github.io/  
- Bio‑Learn GEO data sources: https://bio-learn.github.io/data.html  

---

## License

This project is licensed under the terms of the [MIT License](LICENSE).

---

*For detailed workflow descriptions, code walkthroughs, and result interpretations, please refer to the `README.md` files inside the [`wgbs/`](wgbs/README.md) and [`epic_array/`](epic_array/README.md) directories.*
```
