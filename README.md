# RNA-seq-analysis

A reproducible, Mac-native pipeline to quantify differential gene-expression in BT-549 and MDA-MB-231 breast-cancer cell lines after **BRCA1 knock-down (KD)** and/or **TNFα stimulation**.

---

## ✨ Key Features

| Stage                     | Tool                     | Notes                                          |
| ------------------------- | ------------------------ | ---------------------------------------------- |
| QC                        | **FastQC & MultiQC**     | Aggregated HTML report for all 48 FASTQ files  |
| Transcript quantification | **Kallisto**             | 100 bootstraps; skips if results already exist |
| Differential expression   | **DESeq2 (R)**           | Contrast = `KD` vs `CTRL` (modifiable)         |
| Visualisation             | **Matplotlib / Seaborn** | PCA, volcano plot, heat-map (top genes)        |
| Annotation                | **Ensembl GRCh38 cDNA**  | Index built once, reused for all samples       |

---

## 🔖 Directory Layout

```
project/
├── analysis.py              # main pipeline
├── fastqc/                  # raw data: one dir per sample (FASTQ.gz)
├── refs/                    # reference files & kallisto index
├── quant/                   # Kallisto outputs (auto-created)
├── results/                 # counts + DE tables (auto-created)
└── plots/                   # figures (auto-created)
```

---

## ⚙️ Requirements

| Category     | Version / Package                              |
| ------------ | ---------------------------------------------- |
| Python ≥ 3.8 | `matplotlib pandas numpy seaborn scikit-learn` |
| R runtime    | `Rscript`, **DESeq2 ≥ 1.38**, `optparse`       |
| CLI tools    | `fastqc`, `multiqc`, `kallisto`                |

<details>
<summary>Quick install (Homebrew + pip)</summary>

```bash
brew install fastqc kallisto r
pip install -U matplotlib pandas numpy seaborn scikit-learn
Rscript -e 'if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org/"); BiocManager::install("DESeq2")'
Rscript -e 'install.packages("optparse", repos="https://cloud.r-project.org/")'
```

</details>

---

## 🚀 Quick Start

```bash
# 1 — clone / copy project
cd project

# 2 — download human cDNA & build Kallisto index (one-off)
cd refs
curl -O ftp://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
kallisto index -i grch38-cdna-kallisto.idx Homo_sapiens.GRCh38.cdna.all.fa.gz
cd ..

# 3 — place FASTQ files → fastqc/SAMPLE/
#     naming examples: CTRL-BT-BSA-1, sh-MDA-TNF-24

# 4 — run pipeline
python3 analysis.py
```

Outputs:

* `results/` → counts (`gene_counts.tsv`), DE table (`DE_KD_vs_CTRL.tsv`), MultiQC report
* `plots/`   → `pca.png`, `volcano_kd_vs_ctrl.png`, `heatmap_top50.png`

---

## 🔄 Re-running / Partial Runs

* **Kallisto**: skips if `quant/SAMPLE/abundance.tsv` exists.
* **FastQC**: overwrite old QC by deleting `fastqc/*/qc/` before re-running.

---

## 📝 Citation

If you use this pipeline, please cite:

* Bray *et al.*, *Nature Biotech.* **34**, 525–527 (2016) — Kallisto
* Love *et al.*, *Genome Biol.* **15**, 550 (2014) — DESeq2

---

© 2025 Sean Posner