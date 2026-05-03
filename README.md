# 🧬 Epigenetics and Aging Analysis

This repository contains two fully documented analyses exploring DNA methylation biology:

- **Whole-Genome Bisulfite Sequencing (WGBS)** pipeline on breast cancer methylomes (based on Lin *et al.*, 2015)
- **Epigenetic Aging Clock Benchmarking** using the Bio-Learn library across public 450K/EPIC datasets

Both analyses include detailed biological context, tool explanations, reproducible workflows, and interpretation of results.

> 📖 See individual subdirectory READMEs for full walkthroughs and explanations.

---

## 👤 Author

**Asim Ahmed**  
BSBI-2023 — CMS ID: 454572  
NUST SINES  

**Due:** 3 May 2026

---

## 📁 Repository Structure

```
.
├── wgbs/
├── epic_array/
├── README.md
├── LICENSE
├── requirements.txt
└── .gitignore
```

---

## ⚡ Quick Start

```bash
pip install -r requirements.txt
USE_REAL_DATA=1 python epic_array/src/analysis.py
```

---

## 📜 License

MIT
