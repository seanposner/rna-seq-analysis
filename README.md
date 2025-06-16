# RNA-seq-analysis

A reproducible, Mac-native pipeline for quantifying differential gene expression in BT-549 and MDA-MB-231 breast cancer cell lines after **BRCA1 knockdown (KD)** and/or **TNFα stimulation**.

---

## ✨ Key Features

| Stage                   | Tool/Resource       | Notes                                      |
| ----------------------- | ------------------- | ------------------------------------------ |
| QC                      | FastQC, MultiQC     | Aggregated HTML report for all FASTQ files |
| UMI Extraction          | UMI-tools           | UMI-aware deduplication (if UMIs present)  |
| Alignment               | STAR                | Fast, splice-aware alignment to GRCh38     |
| Counting                | featureCounts       | Gene-level count matrix                    |
| Differential Expression | DESeq2 (R)          | Customizable contrasts (e.g., KD vs CTRL)  |
| Visualization           | Matplotlib, Seaborn | PCA, volcano, and heatmap (top genes)      |
| Annotation              | GENCODE/Ensembl     | Human genome reference & annotation        |

---

## 📁 Directory Structure

```
project/
├── analysis.py                 # Main pipeline
├── fastqc/                     # One dir per sample (FASTQ.gz)
├── refs/                       # Reference data & indices
│   ├── grch38-cdna-kallisto.idx
│   ├── star_index/             # STAR genome index (prebuilt)
│   └── gencode.v36.annotation.gtf
├── quant/                      # Per-sample quantification output
├── results/                    # Counts, DE tables
└── plots/                      # Figures (auto-created)
```

---

## ⚙️ Requirements

| Category    | Package/Resource                                                              |
| ----------- | ----------------------------------------------------------------------------- |
| Python ≥3.8 | matplotlib, pandas, numpy, seaborn, scikit-learn                              |
| R           | Rscript, DESeq2 (≥1.38), optparse                                             |
| CLI tools   | fastqc, multiqc, kallisto, STAR, umi-tools, samtools, subread (featureCounts) |

### Installation (macOS/Homebrew)

```bash
brew install fastqc kallisto star umi-tools samtools subread r
pip install -U matplotlib pandas numpy seaborn scikit-learn
Rscript -e 'if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org/"); BiocManager::install("DESeq2")'
Rscript -e 'install.packages("optparse", repos="https://cloud.r-project.org/")'
```

---

## 📥 Reference Data & STAR Genome Index

**You may use a prebuilt STAR genome index (recommended):**

### Option 1: Download a Prebuilt Index

```bash
cd refs/

# ENCODE prebuilt STAR index (GRCh38, GENCODE v36; reliable, fast)
wget https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.tar.gz -O grch38_star_index.tar.gz
tar -xzvf grch38_star_index.tar.gz
mv starIndex star_index

# Download annotation
wget https://gdc.cancer.gov/files/public/file/gdc/GRCh38.d1.vd1/gencode.v36.annotation.gtf.gz
gunzip gencode.v36.annotation.gtf.gz
```

### Option 2: Download from Google Cloud Public Datasets

```bash
pip install gsutil  # if not present
gsutil -m cp -r gs://gcp-public-data--broad-references/hg38/v0/star/<STAR_INDEX> ./star_index
```

Where the STAR_INDEX can be:
```
star_2.7.9a_primary_gencode_human_v27.tar
```

### Option 3: Build Your Own (Advanced)

If you wish to build the STAR index yourself:

```bash
STAR --runThreadN 4 \
     --runMode genomeGenerate \
     --genomeDir star_index \
     --genomeFastaFiles GRCh38.primary_assembly.genome.fa \
     --sjdbGTFfile gencode.v36.annotation.gtf \
     --sjdbOverhang 99
```

---

## 🚀 Quick Start

```bash
# 1 — Clone/copy the project
cd project

# 2 — Prepare references as above (refs/)

# 3 — Place your FASTQ files → fastqc/SAMPLE/
#     Example: fastqc/CTRL-BT-BSA-1/CTRL-BT-BSA-1_R1.fastq.gz

# 4 — Run the pipeline
python3 analysis.py
```

Outputs:

* `results/` → gene\_counts.tsv, DE table, MultiQC report
* `plots/`   → pca.png, volcano\_kd\_vs\_ctrl.png, heatmap\_top50.png

---

## 🔄 Re-running / Partial Runs

* **STAR + UMI-tools**: skips quantification if output exists in `quant/SAMPLE/SAMPLE_counts.txt`
* **FastQC/MultiQC**: QC step is skipped if `results/multiqc/multiqc_report.html` exists
* **Overwrite QC**: Delete `fastqc/*/qc/` to rerun QC

---

## 📝 Citation

If you use this pipeline, please cite:

* Bray *et al.*, Nat. Biotech. **34**, 525–527 (2016) — Kallisto
* Dobin *et al.*, Bioinformatics **29**, 15–21 (2013) — STAR
* Smith *et al.*, Genome Res. **27**, 491–499 (2017) — UMI-tools
* Liao *et al.*, Bioinformatics **30**, 923–930 (2014) — featureCounts
* Love *et al.*, Genome Biol. **15**, 550 (2014) — DESeq2

---

© 2025 Sean Posner