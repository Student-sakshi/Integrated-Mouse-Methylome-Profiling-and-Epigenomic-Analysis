#!/bin/bash

# ============================================================
# Mouse Methylome Analysis
# WGBS preprocessing and methylation calling
#
# Samples:
#   Mouse kidney and liver
#   Developmental stages E14.5 and E15.5
#   2 biological replicates per tissue/timepoint
#
# Reference:
#   mm10 / GRCm38 reduced reference (chr18-chr19)
#
# Workflow:
#   FASTQ
#   -> FastQC
#   -> Trim Galore!
#   -> Bismark alignment
#   -> Deduplication
#   -> Methylation extraction
#   -> Qualimap / Preseq
#   -> MultiQC
# ============================================================

set -euo pipefail

# ------------------------------------------------------------
# Project paths
# ------------------------------------------------------------

PROJECT="/vol/COMPEPIWS/groups/wgbs3"
TASKS="${PROJECT}/tasks"
DATA="${PROJECT}/data"
RESULTS="${TASKS}/results"

REFERENCE="/vol/COMPEPIWS/pipelines/references"
CONFIG="/vol/COMPEPIWS/pipelines/config"

mkdir -p "${TASKS}"
mkdir -p "${DATA}"
mkdir -p "${RESULTS}"

# ------------------------------------------------------------
# Input data
# ------------------------------------------------------------

RAW_DATA="/vol/COMPEPIWS/data/reduced/WGBS"

# The dataset contains eight single-end WGBS samples:
#
# kidney_14.5.1
# kidney_14.5.2
# kidney_15.5.1
# kidney_15.5.2
# liver_14.5.1
# liver_14.5.2
# liver_15.5.1
# liver_15.5.2

# Create symbolic links to the raw FASTQ files
cd "${DATA}"

ln -sf "${RAW_DATA}/kidney_14.5_WGBS_1_1.fastq.gz" .
ln -sf "${RAW_DATA}/kidney_14.5_WGBS_2_1.fastq.gz" .
ln -sf "${RAW_DATA}/kidney_15.5_WGBS_1_1.fastq.gz" .
ln -sf "${RAW_DATA}/kidney_15.5_WGBS_2_1.fastq.gz" .
ln -sf "${RAW_DATA}/liver_14.5_WGBS_1_1.fastq.gz" .
ln -sf "${RAW_DATA}/liver_14.5_WGBS_2_1.fastq.gz" .
ln -sf "${RAW_DATA}/liver_15.5_WGBS_1_1.fastq.gz" .
ln -sf "${RAW_DATA}/liver_15.5_WGBS_2_1.fastq.gz"

# ------------------------------------------------------------
# Sample sheet
# ------------------------------------------------------------

SAMPLESHEET="${TASKS}/samplesheet.csv"

cat > "${SAMPLESHEET}" <<EOF
sample,fastq_1,fastq_2
kidney_14.5.1,${DATA}/kidney_14.5_WGBS_1_1.fastq.gz,
kidney_14.5.2,${DATA}/kidney_14.5_WGBS_2_1.fastq.gz,
kidney_15.5.1,${DATA}/kidney_15.5_WGBS_1_1.fastq.gz,
kidney_15.5.2,${DATA}/kidney_15.5_WGBS_2_1.fastq.gz,
liver_14.5.1,${DATA}/liver_14.5_WGBS_1_1.fastq.gz,
liver_14.5.2,${DATA}/liver_14.5_WGBS_2_1.fastq.gz,
liver_15.5.1,${DATA}/liver_15.5_WGBS_1_1.fastq.gz,
liver_15.5.2,${DATA}/liver_15.5_WGBS_2_1.fastq.gz,
EOF

# ------------------------------------------------------------
# Load computational environment
# ------------------------------------------------------------

source /vol/COMPEPIWS/conda/miniconda3/bin/activate \
/vol/COMPEPIWS/conda/miniconda3/envs/core

# ------------------------------------------------------------
# Run nf-core/methylseq
# ------------------------------------------------------------

nextflow run nf-core/methylseq \
    -r 2.3.0 \
    -profile singularity \
    -params-file "${CONFIG}/wgbs_params.json" \
    --input "${SAMPLESHEET}" \
    --outdir "${RESULTS}" \
    -process.maxForks 2

# ------------------------------------------------------------
# Main processing performed by the workflow
#
# 1. Reference genome preparation
# 2. Raw-read quality control
# 3. Adapter trimming with Trim Galore!
# 4. Bisulfite-aware alignment using Bismark/Bowtie2
# 5. Duplicate removal
# 6. Methylation extraction
# 7. Bismark reporting
# 8. Qualimap alignment QC
# 9. Preseq complexity estimation
# 10. MultiQC summary
# ------------------------------------------------------------

echo "WGBS preprocessing completed."
echo "Results: ${RESULTS}"
