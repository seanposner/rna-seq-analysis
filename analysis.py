"""
RNA-seq Analysis Pipeline: BRCA1 Knockdown ± TNFalpha in BT549 & MDA-MB-231 Cells (Mac Native)
==========================================================================================

Overview:
    - End-to-end RNA-seq workflow for quantifying transcriptome changes across 8 biological conditions.
    - Handles QC, quantification, DE analysis (via R/DESeq2), and result visualization.

Pipeline Stages:
    1. **Quality Control (QC):**
        - Runs FastQC on each FASTQ file.
        - Aggregates QC metrics using MultiQC.
    2. **Quantification:**
        - Runs Kallisto (paired-end, transcript-level) for each sample.
        - Skips quantification if output already exists.
    3. **Collation:**
        - Extracts counts and TPM matrices from Kallisto output.
        - Compiles per-sample metadata (parsed from folder names).
    4. **Differential Expression:**
        - Runs DESeq2 via Rscript (requires R/DESeq2/optparse).
        - Outputs a table of log2 fold changes, p-values, etc.
    5. **Visualization:**
        - Generates PCA, volcano plots, and clustered heatmaps for key results.

Folder Structure Example:
    Analysis/
        analysis.py
        fastqc/
            CTRL-BT-BSA-1/    # (Contains .fastq.gz files)
            sh-MDA-TNF-22/
            ...
        refs/
            grch38-cdna-kallisto.idx
        results/
        plots/
        quant/

Author: Sean Posner
License: MIT
"""

import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA

ROOT          = Path(__file__).resolve().parent
THREADS       = int(os.getenv("THREADS", os.cpu_count() or 4))
KALLISTO_IDX  = ROOT / "refs" / "grch38-cdna-kallisto.idx"
RESULTS_DIR   = ROOT / "results"
PLOTS_DIR     = ROOT / "plots"
FASTQC_DIR    = ROOT / "fastqc"
QUANT_DIR     = ROOT / "quant"

FASTQC_EXE   = shutil.which("fastqc")   or "fastqc"
MULTIQC_EXE  = shutil.which("multiqc")  or "multiqc"
KALLISTO_EXE = shutil.which("kallisto") or "kallisto"

# Updated regex to handle both "CTRL-BT-BSA-1" and "sh-BT-BSA-7"
SAMPLE_ID_PATTERN = re.compile(
    r"^(CTRL|KD|sh)[-_](BT|MDA)[-_](BSA|TNF)[-_](\d+)$", re.IGNORECASE
)

def _assert_tool(name: str, exe: str):
    if shutil.which(exe) is None:
        sys.exit(f"[ERROR] Required tool '{name}' not found in $PATH.")

for tool, exe in {
    "FastQC": FASTQC_EXE,
    "MultiQC": MULTIQC_EXE,
    "Kallisto": KALLISTO_EXE,
    "Rscript": shutil.which("Rscript") or "Rscript"
}.items():
    _assert_tool(tool, exe)

def _run(cmd: List[str]):
    try:
        subprocess.run([str(x) for x in cmd], check=True)
    except subprocess.CalledProcessError as e:
        sys.exit(f"\n[FAIL] Command failed ({e.returncode}): {' '.join(map(str, e.cmd))}\n")

def fastqc_dirwise(sample_dir: Path, outdir: Path):
    fqs = sorted(sample_dir.glob("*.fastq.gz"))
    if not fqs:
        print(f"[SKIP] No FASTQ files in {sample_dir}")
        return
    outdir.mkdir(parents=True, exist_ok=True)
    _run([FASTQC_EXE, "--quiet", "--threads", str(THREADS), "--outdir", str(outdir)] + [str(f) for f in fqs])

def multiqc_on_qc(qc_dir: Path, out: Path):
    out.mkdir(parents=True, exist_ok=True)
    _run([MULTIQC_EXE, str(qc_dir), "--outdir", str(out)])

def kallisto_quant(sample: str, fq_dir: Path, out: Path):
    out.mkdir(parents=True, exist_ok=True)
    fq_files = sorted(fq_dir.glob("*.fastq.gz"))
    # Only run Kallisto if not already run (abundance.tsv exists)
    abun_file = out / "abundance.tsv"
    if abun_file.exists():
        print(f"[SKIP] Quant already exists for {sample}")
        return
    if len(fq_files) == 2:
        cmd = [
            str(KALLISTO_EXE), "quant",
            "-i", str(KALLISTO_IDX),
            "-o", str(out),
            "-b", "100",
            "-t", str(THREADS),
            str(fq_files[0]), str(fq_files[1])
        ]
    else:
        sys.exit(f"[ERROR] Expected 2 FASTQ files for paired-end: {fq_dir}")
    _run(cmd)

def parse_sample_id(sample_dir_name: str) -> Dict[str, str]:
    """
    Parse sample folder name into metadata.
    Handles both "CTRL-BT-BSA-1" and "sh-BT-BSA-7" formats.
    """
    m = SAMPLE_ID_PATTERN.match(sample_dir_name)
    if not m:
        raise ValueError(
            f"Bad sample ID '{sample_dir_name}'. Expected: "
            "CTRL|KD|sh-BT|MDA-BSA|TNF-#"
        )
    return {
        "sample": sample_dir_name,
        "brca1_kd": m.group(1).upper(),      # CTRL, KD, or sh
        "cell_line": m.group(2).upper(),     # BT, MDA
        "tnf": m.group(3).upper(),           # BSA, TNF
        "rep": m.group(4)
    }

def import_kallisto_quant(qdir: Path) -> pd.DataFrame:
    df = pd.read_csv(qdir / "abundance.tsv", sep="\t").set_index("target_id")
    return df[["est_counts", "tpm"]].rename(columns={"est_counts": "counts"})

def build_matrices(quant_dir: Path):
    counts, tpms, meta_rows = [], [], []
    detected = []
    for qdir in sorted(quant_dir.iterdir()):
        abun = qdir / "abundance.tsv"
        if not abun.exists():
            print(f"[WARN] Skipping {qdir.name} (missing abundance.tsv)")
            continue
        sample = qdir.name
        try:
            meta = parse_sample_id(sample)
            detected.append(sample)
        except Exception as e:
            print(f"[SKIP] {sample}: {e}")
            continue
        meta_rows.append(meta)
        try:
            df = import_kallisto_quant(qdir)
            counts.append(df["counts"].rename(sample))
            tpms.append(df["tpm"].rename(sample))
        except Exception as e:
            print(f"[WARN] Failed to load quant for {sample}: {e}")
    print(f"\n[INFO] {len(detected)} samples with quantification found:")
    for s in detected:
        print(f"  - {s}")
    if not counts:
        sys.exit("[ERROR] No valid Kallisto quantification results found in quant/. Aborting.")
    meta = pd.DataFrame(meta_rows).set_index("sample")
    return pd.concat(counts, axis=1), pd.concat(tpms, axis=1), meta

# ---- Differential Expression with R/DESeq2 ----

DESEQ2_R_SCRIPT = ROOT / "deseq2_deseq_contrast.R"

def write_deseq2_r():
    # Write an R script to run DESeq2 DE analysis
    DESEQ2_R_SCRIPT.write_text("""
suppressMessages(library("optparse"))
suppressMessages(library("DESeq2"))

option_list <- list(
  make_option(c("--counts"), type="character", help="counts matrix tsv"),
  make_option(c("--meta"), type="character", help="sample meta tsv"),
  make_option(c("--factor"), type="character", help="test factor"),
  make_option(c("--test"), type="character", help="test level"),
  make_option(c("--ref"), type="character", help="ref level"),
  make_option(c("--out"), type="character", help="output tsv")
)
opt <- parse_args(OptionParser(option_list=option_list))

cts <- as.matrix(read.table(opt$counts, sep="\\t", header=TRUE, row.names=1, check.names=FALSE))
meta <- read.table(opt$meta, sep="\\t", header=TRUE, row.names=1, check.names=FALSE)
meta[] <- lapply(meta, as.factor)

dds <- DESeqDataSetFromMatrix(countData=cts, colData=meta, design=as.formula(paste("~ cell_line + tnf +", opt$factor)))
dds[[opt$factor]] <- relevel(dds[[opt$factor]], ref=opt$ref)
dds <- DESeq(dds)
res <- results(dds, contrast=c(opt$factor, opt$test, opt$ref))
res <- as.data.frame(res)
res$gene <- rownames(res)
write.table(res, file=opt$out, sep="\\t", quote=FALSE, row.names=FALSE)
    """.strip())

def run_deseq2_r(counts_path, meta_path, factor, test, ref, out_path):
    write_deseq2_r()
    cmd = [
        "Rscript", str(DESEQ2_R_SCRIPT),
        "--counts", str(counts_path),
        "--meta", str(meta_path),
        "--factor", factor,
        "--test", test,
        "--ref", ref,
        "--out", str(out_path)
    ]
    _run(cmd)

def pca_plot(vst: pd.DataFrame, meta: pd.DataFrame, out: Path):
    from sklearn.decomposition import PCA
    pcs: np.ndarray = PCA(n_components=2).fit_transform(vst.T)
    fig, ax = plt.subplots(figsize=(5, 4))
    for (brca, tnf, line), grp in meta.groupby(["brca1_kd", "tnf", "cell_line"]):
        idx = np.array([vst.columns.get_loc(s) for s in grp.index])
        ax.scatter(pcs[idx, 0], pcs[idx, 1], label=f"{line}-{brca}-{tnf}", s=60, alpha=.75)
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.legend(fontsize="x-small", bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)

def volcano_plot(de: pd.DataFrame, out: Path, alpha: float = 0.05):
    de = de.copy()
    de["nlp10"] = -np.log10(de["padj"].clip(1e-300))
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(de["log2FoldChange"], de["nlp10"], s=8, alpha=.5)
    ax.axhline(-np.log10(alpha), linestyle="--", color="grey")
    ax.set_xlabel("log2 fold‑change")
    ax.set_ylabel("-log10 adj‑p")
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)

def heatmap(vst: pd.DataFrame, meta: pd.DataFrame, genes: List[str], out: Path):
    data = vst.loc[genes].dropna()
    annot = meta[["cell_line", "brca1_kd", "tnf"]].copy()
    lut = {'CTRL': 'blue', 'KD': 'red', 'BSA': 'grey', 'TNF': 'orange', 'BT': 'green', 'MDA': 'purple'}

    # Apply color mapping column-wise and ensure result is a DataFrame
    col_colors = annot.apply(lambda col: col.map(lambda x: lut.get(x, "black")))

    g = sns.clustermap(data, cmap="vlag", z_score=0, col_colors=col_colors)
    plt.savefig(out, dpi=300)
    plt.close()

def vst_log2(counts: pd.DataFrame) -> pd.DataFrame:
    """Variance-stabilizing log2 transform, DataFrame-preserving."""
    return pd.DataFrame(np.log2(counts + 1), index=counts.index, columns=counts.columns)

def main():
    RESULTS_DIR.mkdir(exist_ok=True)
    PLOTS_DIR.mkdir(exist_ok=True)
    QUANT_DIR.mkdir(exist_ok=True)

    print("[STEP] QC")
    multiqc_report = RESULTS_DIR / "multiqc" / "multiqc_report.html"
    if multiqc_report.exists():
        print(f"[SKIP] FastQC & MultiQC already completed: {multiqc_report}")
    else:
        for sample_dir in sorted(FASTQC_DIR.iterdir()):
            if not sample_dir.is_dir():
                continue
            fastqc_dirwise(sample_dir, sample_dir / "qc")
        multiqc_on_qc(FASTQC_DIR, RESULTS_DIR / "multiqc")

    print("[STEP] Kallisto quant")
    for sample_dir in sorted(FASTQC_DIR.iterdir()):
        if not sample_dir.is_dir():
            continue
        outdir = QUANT_DIR / sample_dir.name
        kallisto_quant(sample_dir.name, sample_dir, outdir)

    print("[STEP] Collate counts")
    counts, tpms, meta = build_matrices(QUANT_DIR)
    counts = counts.round().astype(int)
    counts_tsv = RESULTS_DIR / "gene_counts.tsv"
    tpms_tsv = RESULTS_DIR / "gene_tpms.tsv"
    meta_tsv = RESULTS_DIR / "sample_meta.tsv"
    counts.to_csv(counts_tsv, sep="\t")
    tpms.to_csv(tpms_tsv, sep="\t")
    meta.to_csv(meta_tsv, sep="\t")

    print("[STEP] Differential expression")
    de_out = RESULTS_DIR / "DE_KD_vs_CTRL.tsv"
    run_deseq2_r(counts_tsv, meta_tsv, "brca1_kd", "KD", "CTRL", de_out)
    de_kd = pd.read_csv(de_out, sep="\t")

    vst = vst_log2(counts)
    print("[STEP] Figures")
    pca_plot(vst, meta, PLOTS_DIR / "pca.png")
    volcano_plot(de_kd, PLOTS_DIR / "volcano_kd_vs_ctrl.png")
    sig_genes = de_kd.query("padj < 0.05 & abs(log2FoldChange) > 1")["gene"].tolist()[:50]
    if sig_genes:
        heatmap(vst, meta, sig_genes, PLOTS_DIR / "heatmap_top50.png")

    print("[DONE] Results in ./results and ./plots ✓")

if __name__ == "__main__":
    main()