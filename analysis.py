"""
RNA-seq Analysis Pipeline: BRCA1 Knockdown ± TNFalpha in BT549 & MDA-MB-231 Cells (Mac Native)
==========================================================================================

Overview:
    - End-to-end RNA-seq workflow for quantifying transcriptome changes across 8 biological conditions.
    - Handles QC, quantification, DE analysis (via R/DESeq2), and result visualization.

Pipeline Stages:
    1. Quality Control (QC): FastQC + MultiQC (skipped if already done)
    2. Quantification: Kallisto (paired-end, transcript-level, skipped if abundance.tsv exists)
    3. Collation: Extract count/TPM matrices and per-sample metadata (from folder names)
    4. Differential Expression: DESeq2 in R (requires optparse, DESeq2)
    5. Visualization: PCA, volcano, heatmaps

Folder Structure Example:
    Analysis/
        analysis.py
        fastqc/
            CTRL-BT-BSA-1/    # Contains .fastq.gz files
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
from typing import Dict, List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA

ROOT          = Path(__file__).resolve().parent
THREADS       = 4
KALLISTO_IDX  = ROOT / "refs" / "grch38-cdna-kallisto.idx"
RESULTS_DIR   = ROOT / "results"
PLOTS_DIR     = ROOT / "plots"
FASTQC_DIR    = ROOT / "fastqc"
QUANT_DIR     = ROOT / "quant"

FASTQC_EXE   = shutil.which("fastqc")   or "fastqc"
MULTIQC_EXE  = shutil.which("multiqc")  or "multiqc"
KALLISTO_EXE = shutil.which("kallisto") or "kallisto"

# Accepts both CTRL-BT-BSA-1 and sh-BT-BSA-7 and KD-BT-BSA-7
SAMPLE_ID_PATTERN = re.compile(
    r"^(CTRL|KD|SH)[-_](BT|MDA)[-_](BSA|TNF)[-_](\d+)$", re.IGNORECASE
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

def umi_extract_and_dedup(sample_name, r1_path, r2_path, star_index, output_dir, annotation_gtf, threads=4):
    """
    UMI-aware processing of a sample.
    - Extracts UMIs using UMI-tools from R2, appends to R1 header.
    - Aligns reads with STAR.
    - Deduplicates reads using UMI-tools.
    - Counts features using featureCounts.
    """
    from pathlib import Path
    import subprocess

    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)
    extracted_r1 = output_dir / f"{sample_name}_R1_extracted.fastq.gz"
    aligned_bam = output_dir / f"{sample_name}_Aligned.sortedByCoord.out.bam"
    dedup_bam = output_dir / f"{sample_name}_dedup.bam"
    counts_txt = output_dir / f"{sample_name}_counts.txt"

    # 1. UMI extraction
    subprocess.run([
        "umi_tools", "extract",
        "--bc-pattern=NNNNNNNN",  # 8bp UMI, update if needed
        "-I", str(r1_path),
        "--read2-in", str(r2_path),
        "-S", str(extracted_r1)
    ], check=True)

    # 2. STAR alignment
    subprocess.run([
        "STAR",
        "--runThreadN", str(threads),
        "--genomeDir", str(star_index),
        "--readFilesIn", str(extracted_r1),
        "--readFilesCommand", "gunzip", "-c",
        "--outFileNamePrefix", str(output_dir / f"{sample_name}_"),
        "--outSAMtype", "BAM", "SortedByCoordinate"
    ], check=True)

    star_bam = output_dir / f"{sample_name}_Aligned.sortedByCoord.out.bam"
    
    # 2.5. Index the BAM
    subprocess.run([
        "samtools", "index", str(star_bam)
    ], check=True)

    # 3. Deduplication with UMI-tools
    subprocess.run([
        "umi_tools", "dedup",
        "-I", str(star_bam),
        "-S", str(dedup_bam)
    ], check=True)

    # 4. Quantification with featureCounts
    subprocess.run([
        "featureCounts",
        "-T", str(threads),
        "-a", str(annotation_gtf),
        "-o", str(counts_txt),
        str(dedup_bam)
    ], check=True)

    return counts_txt

def parse_sample_id(sample_dir_name: str) -> Dict[str, str]:
    """
    Parse sample folder name into metadata.
    Handles both "CTRL-BT-BSA-1", "KD-BT-BSA-7", and "sh-BT-BSA-7" formats.
    Always returns 'KD' for KD/SH, and 'CTRL' for controls.
    """
    m = SAMPLE_ID_PATTERN.match(sample_dir_name)
    if not m:
        raise ValueError(
            f"Bad sample ID '{sample_dir_name}'. Expected: "
            "CTRL|KD|SH-BT|MDA-BSA|TNF-#"
        )
    group = m.group(1).upper()
    # Map 'SH' to 'KD' for analysis consistency
    brca1_kd = "KD" if group in ("KD", "SH") else "CTRL"
    return {
        "sample": sample_dir_name,
        "brca1_kd": brca1_kd,                # Only "CTRL" or "KD"
        "cell_line": m.group(2).upper(),     # BT, MDA
        "tnf": m.group(3).upper(),           # BSA, TNF
        "rep": m.group(4)
    }

# ---- Differential Expression with R/DESeq2 ----

DESEQ2_R_SCRIPT = ROOT / "deseq2_deseq_contrast.R"

def write_deseq2_r():
    DESEQ2_R_SCRIPT.write_text("""
suppressMessages(library("optparse"))
suppressMessages(library("DESeq2"))

option_list <- list(
  make_option(c("--counts"), type="character", help="counts matrix tsv"),
  make_option(c("--meta"), type="character", help="sample meta tsv"),
  make_option(c("--factor"), type="character", help="test factor"),
  make_option(c("--test"), type="character", help="test level"),
  make_option(c("--ref"), type="character", help="ref level"),
  make_option(c("--deout"), type="character", help="output DE tsv"),
  make_option(c("--vstout"), type="character", help="output VST tsv")
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
write.table(res, file=opt$deout, sep="\\t", quote=FALSE, row.names=FALSE)

# VST matrix
vst <- vst(dds)
vst_mat <- assay(vst)
write.table(vst_mat, file=opt$vstout, sep="\\t", quote=FALSE)
    """.strip())

def run_deseq2_r(counts_path, meta_path, factor, test, ref, de_out_path, vst_out_path):
    write_deseq2_r()
    cmd = [
        "Rscript", str(DESEQ2_R_SCRIPT),
        "--counts", str(counts_path),
        "--meta", str(meta_path),
        "--factor", factor,
        "--test", test,
        "--ref", ref,
        "--deout", str(de_out_path),
        "--vstout", str(vst_out_path)
    ]
    _run(cmd)

def pca_plot(vst: pd.DataFrame, meta: pd.DataFrame, out: Path):
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
    ax.set_xlabel("log2 fold-change")
    ax.set_ylabel("-log10 adj-p")
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)

def heatmap(vst: pd.DataFrame, meta: pd.DataFrame, genes: List[str], out: Path):
    data = vst.loc[genes].dropna()
    annot = meta[["cell_line", "brca1_kd", "tnf"]].copy()
    lut = {'CTRL': 'blue', 'KD': 'red', 'BSA': 'grey', 'TNF': 'orange', 'BT': 'green', 'MDA': 'purple'}
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

    # ---- Paths to STAR and GTF ----
    STAR_INDEX = ROOT / "refs" / "star"
    ANNOTATION_GTF = ROOT / "refs" / "gencode.v45.annotation.gtf"

    print("[STEP] QC")
    multiqc_report = RESULTS_DIR / "multiqc" / "multiqc_report.html"
    if multiqc_report.exists():
        print(f"[SKIP] FastQC & MultiQC already completed: {multiqc_report})")
    else:
        for sample_dir in sorted(FASTQC_DIR.iterdir()):
            if not sample_dir.is_dir():
                continue
            fastqc_dirwise(sample_dir, sample_dir / "qc")
        multiqc_on_qc(FASTQC_DIR, RESULTS_DIR / "multiqc")

    print("[STEP] UMI-aware quantification (STAR + UMI-tools + featureCounts)")
    counts_files = []
    meta_rows = []
    for sample_dir in sorted(FASTQC_DIR.iterdir()):
        if not sample_dir.is_dir():
            continue
        sample_name = sample_dir.name
        outdir = QUANT_DIR / sample_name
        counts_txt = outdir / f"{sample_name}_counts.txt"
        if counts_txt.exists():
            print(f"[SKIP] {sample_name}: already quantified.")
            counts_files.append((sample_name, counts_txt))
            try:
                meta = parse_sample_id(sample_name)
                meta_rows.append(meta)
            except Exception as e:
                print(f"[SKIP META] {sample_name}: {e}")
            continue
        fq_r1 = sorted(sample_dir.glob("*R1*.fastq.gz"))
        fq_r2 = sorted(sample_dir.glob("*R2*.fastq.gz"))
        if len(fq_r1) != 1 or len(fq_r2) != 1:
            print(f"[WARN] Skipping {sample_name}: Need exactly 1 R1 and 1 R2 FASTQ")
            continue
        counts_txt = umi_extract_and_dedup(
            sample_name, fq_r1[0], fq_r2[0],
            STAR_INDEX, outdir, ANNOTATION_GTF, threads=THREADS
        )
        counts_files.append((sample_name, counts_txt))
        try:
            meta = parse_sample_id(sample_name)
            meta_rows.append(meta)
        except Exception as e:
            print(f"[SKIP META] {sample_name}: {e}")

    print("[STEP] Collate counts")
    all_counts = []
    for sample_name, counts_txt in counts_files:
        # Parse featureCounts output (skip header lines, grab first sample column)
        df = pd.read_csv(counts_txt, sep='\t', comment='#')
        if sample_name not in df.columns:
            # fallback: try last column (featureCounts sometimes labels with full path)
            sample_col = df.columns[-1]
        else:
            sample_col = sample_name
        df = df.set_index('Geneid')[[sample_col]]
        df.columns = [sample_name]
        all_counts.append(df)
    if not all_counts:
        sys.exit("[ERROR] No quantification results found. Aborting.")
    counts = pd.concat(all_counts, axis=1)
    counts = counts.loc[~counts.index.str.startswith('__')]  # Remove summary rows

    meta = pd.DataFrame(meta_rows).set_index("sample")
    counts_tsv = RESULTS_DIR / "gene_counts.tsv"
    meta_tsv = RESULTS_DIR / "sample_meta.tsv"
    counts.to_csv(counts_tsv, sep="\t")
    meta.to_csv(meta_tsv, sep="\t")

    print("[STEP] Differential expression + VST")
    de_out = RESULTS_DIR / "DE_KD_vs_CTRL.tsv"
    vst_out = RESULTS_DIR / "vst_counts.tsv"
    run_deseq2_r(counts_tsv, meta_tsv, "brca1_kd", "KD", "CTRL", de_out, vst_out)
    de_kd = pd.read_csv(de_out, sep="\t")
    vst = pd.read_csv(vst_out, sep="\t", index_col=0)

    # Filter to variable genes for PCA (optional but helps!)
    variances = vst.var(axis=1)
    top_genes = variances.nlargest(500).index
    vst_top = vst.loc[top_genes]

    print("[STEP] Figures")
    pca_plot(vst_top, meta, PLOTS_DIR / "pca.png")
    volcano_plot(de_kd, PLOTS_DIR / "volcano_kd_vs_ctrl.png")
    sig_genes = de_kd.query("padj < 0.05 & abs(log2FoldChange) > 1")["gene"].tolist()[:50]
    if sig_genes:
        heatmap(vst, meta, sig_genes, PLOTS_DIR / "heatmap_top50.png")

    print("[DONE] Results in ./results and ./plots ✓")

if __name__ == "__main__":
    main()