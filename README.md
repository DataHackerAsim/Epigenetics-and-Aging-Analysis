🧬 Epigenetics and Aging AnalysisWhole-Genome Bisulfite Sequencing (WGBS) & Epigenetic Aging-Clock BenchmarkingA dual-module bioinformatics pipeline exploring DNA methylation biology through genomic sequencing data and microarray-based aging biomarkers.📌 Project OverviewThis repository contains two fully documented, reproducible analyses covering distinct aspects of DNA methylation biology:Whole-Genome Bisulfite Sequencing (WGBS) Pipeline: Applied to breast-cancer methylomes, following the Galaxy Training Network tutorial and utilizing datasets from Lin et al. (2015).Epigenetic Aging-Clock Benchmarking: Leveraging the Bio-Learn Python library to evaluate 14 distinct aging clocks across two public 450K/EPIC array blood-methylation datasets. The analysis features comprehensive benchmarking via correlation matrices, age-deviation heatmaps, predicted-vs-chronological-age scatters, Mean Absolute Error (MAE) charts, and statistical distribution plots.Author: Asim Ahmed (BSBI-2023, CMS ID: 454572)Institution: NUST SINESTerm: Spring 2026 (Due: 3 May 2026)📑 Table of ContentsRepository StructureModule Breakdown🚀 Quick Start & Reproducibility📊 Methodology & ResultsPart 1: WGBS PipelinePart 2: EPIC Array Benchmarking⚙️ Data Availability & Simulation Mode📚 Citations & References📜 License📂 Repository Structuremethylation_assignment/
│
├── wgbs/                                  # Part 1: WGBS bisulfite sequencing pipeline
│   ├── README.md                          # Detailed workflow + result interpretation
│   ├── run_wgbs_pipeline.sh               # Linux/CLI reproducibility script
│   ├── notes_clustering.md                # Hierarchical-clustering R walkthrough
│   ├── wgbs_workflow.png                  # Pipeline architecture diagram
│   └── key_findings_summary.png           # 3-panel summary of Lin 2015 findings
│
├── epic_array/                            # Part 2: EPIC-array aging-clock benchmarking
│   ├── biolearn_aging_clocks.ipynb        # Primary deliverable notebook (Colab-ready)
│   ├── README.md                          # Detailed workflow + result interpretation
│   ├── src/
│   │   ├── analysis.py                    # Main pipeline (clocks + 5 viz types)
│   │   ├── simulate_data.py               # Methylation simulator (offline fallback)
│   │   └── download_real_data.py          # GSE40279 / GSE41169 fetcher
│   └── results/
│       ├── summary.md                     # Auto-generated metrics report
│       ├── figures/                       # 8 PNGs (3 per-dataset × 2 + 2 cross-dataset)
│       └── tables/                        # Prediction + summary CSVs
│
├── README.md                              # This top-level documentation
├── LICENSE                                # MIT License
├── requirements.txt                       # Python dependencies
└── .gitignore                             # Ignored build/cache files
🧬 Module BreakdownModuleDirectoryPlatform/EnvironmentPrimary Focus1. WGBSwgbs/Galaxy Europe (usegalaxy.eu)Bisulfite QC, sequence alignment, methylation extraction, and DMR detection.2. EPIC Arrayepic_array/Python 3.10+ / Jupyter / biolearn14-clock benchmarking, statistical correlations, heatmaps, and MAE comparison.🚀 Quick Start & ReproducibilityPart 2: EPIC Array (Executable Pipeline)The EPIC Array module is fully executable locally. Ensure you have Python 3.10+ installed.# 1. Clone the repository and navigate to the directory
git clone [https://github.com/yourusername/methylation_assignment.git](https://github.com/yourusername/methylation_assignment.git)
cd methylation_assignment

# 2. Install required dependencies
pip install -r requirements.txt

# 3. Execute the pipeline
# Option A: Run on REAL Bio-Learn datasets (Requires internet; ~500MB download cached to ~/.biolearn/)
USE_REAL_DATA=1 python epic_array/src/analysis.py

# Option B: Run the bundled simulator (For offline or sandboxed environments)
python epic_array/src/analysis.py
Alternatively, you can explore the data interactively via Jupyter Notebook:jupyter lab epic_array/biolearn_aging_clocks.ipynb
Note: Both execution modes automatically regenerate all figures and tables from scratch.Part 1: WGBS (Galaxy Workflows)Part 1 is a documentation-driven deliverable. Because the full WGBS pipeline requires ~500 GB of FASTQ input and a multi-day compute job, it is documented following the Galaxy Training Network format.See wgbs/README.md for the complete walkthrough.See wgbs/run_wgbs_pipeline.sh for the Linux command-line equivalent.📊 Methodology & ResultsPart 1: WGBS PipelineObjective: Apply a whole-genome bisulfite sequencing pipeline to breast-cancer and normal-breast tissue samples from Lin et al. (2015).Architecture: Paired-end bisulfite-converted FASTQ → Falco QC → Trim Galore → bwa-meth → MethylDackel → deepTools profile + Metilene DMR detection → hierarchical HMR clustering across 7 methylomes.Key Findings (Lin et al. summary): (A) Cell lines and tumours demonstrate genome-wide hypomethylation.(B) Bidirectional methylation remodelling occurs (CpG-rich loci hypermethylate, CpG-poor loci hypomethylate).(C) Highly Methylated Regions (HMRs) expand and shift to non-CGI loci in cancer.Part 2: EPIC Array BenchmarkingObjective: Benchmark 14 epigenetic aging clocks across two blood 450K/EPIC array datasets.Clock Correlation: The 13 chronological-age clocks form a tight block ($r > 0.9$), agreeing on relative ordering. As expected, DunedinPACE (a pace-of-aging clock) shows near-zero correlation with absolute chronological-age clocks.Age Deviation: Heatmap sorting samples by chronological age reveals characteristic vertical bands of bias per clock. PhenoAge systematically over-predicts oldest samples (a trait of biological-age clocks).Prediction Accuracy: AltumAge (deep neural net) and Horvath/Hannum classical models achieve excellent accuracy ($r > 0.95$). VidalBralo, StocH, and Horvathv2 achieve the lowest MAE (~10–14 yr).Distribution: Well-calibrated clocks show predicted-age distributions centered cleanly on the dataset's chronological median line.⚙️ Data Availability & Simulation ModeTo ensure maximum reproducibility regardless of execution environment, this repository features a dual-mode execution strategy:Production Mode (USE_REAL_DATA=1): Fetches the real GSE40279 and GSE41169 cohorts directly from NCBI GEO.Simulation Mode (Default): The committed figures in epic_array/results/ were generated using the bundled simulate_data.py module. This was engineered to bypass institutional firewall blocks on ftp.ncbi.nlm.nih.gov.Biological Grounding: The simulator utilizes Bio-Learn's population-mean $\beta$-values and the union of all 14 clocks' coefficient signs to construct $485k \times N$ $\beta$-value matrices. This successfully drives chronological-age clocks to perform with $r > 0.87$ accuracy against ground-truth chronological age, faithfully mimicking real-world microarray behaviors.📚 Citations & ReferencesLiteratureLin, I.-H., Chen, D.-T., Chang, Y.-F., Lee, Y.-L., Su, C.-H., et al. (2015). Hierarchical Clustering of Breast Cancer Methylomes Revealed Differentially Methylated and Expressed Breast Cancer Genes. PLOS ONE 10(2): e0118453. DOI: 10.1371/journal.pone.0118453Ying, K., Paulson, S., Perez-Guevara, M., Emamifar, M., Martinez, M. C., Kwon, D., Poganik, J. R., Moqri, M., Gladyshev, V. N. (2023). Biolearn, an open-source library for biomarkers of aging. bioRxiv. DOI: 10.1101/2023.12.02.569722Wolff, J., Ryan, D., Moosmann, V. (2017). DNA Methylation data analysis. Galaxy Training Materials. LinkTools & DocumentationGalaxy Training Network: DNA Methylation data analysisGalaxy Training Network: Introduction to DNA MethylationBio-Learn Documentation | GEO Data Sources📜 LicenseThis project is licensed under the MIT License. See the LICENSE file for more information.
