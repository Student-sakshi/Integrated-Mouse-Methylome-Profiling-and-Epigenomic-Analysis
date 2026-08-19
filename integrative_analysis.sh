#!/bin/bash

# ============================================================
# Mouse Methylome Analysis
# Integrative epigenomic analysis
#
# Integration:
#   DNA methylation
#   MethylSeekR segmentation
#   ChromHMM
#   ATAC-seq
#   Histone modifications
#   RNA-seq
# ============================================================

set -euo pipefail

BASE="/vol/COMPEPIWS/groups/shared/WGBS/wgbs3"

DIFF="${BASE}/differential"

mkdir -p "${DIFF}"

# ------------------------------------------------------------
# 1. Prepare DMR BED file
# ------------------------------------------------------------

tail -n +2 \
"${DIFF}/mystery_diffmeth_sorted.csv" |
awk -F',' '
BEGIN {
    OFS="\t"
}
{
    gsub(/"/,"",$1)
    print $1,$2,$3
}
' > "${DIFF}/mystery_diffmeth_sorted.bed"

# ------------------------------------------------------------
# 2. Extract heterochromatic ChromHMM states
# ------------------------------------------------------------

grep -E "Het_S|Het_P" \
/vol/COMPEPIWS/groups/shared/ChIP-seq/chipseq1/segmentation/kidney_15_segments.bed \
> "${DIFF}/kidney_Het.bed"

grep -E "Het_S|Het_P" \
/vol/COMPEPIWS/groups/shared/ChIP-seq/chipseq1/segmentation/liver_15_segments.bed \
> "${DIFF}/liver_Het.bed"

# ------------------------------------------------------------
# 3. DMR overlap with heterochromatin
# ------------------------------------------------------------

bedtools intersect \
-a "${DIFF}/mystery_diffmeth_sorted.bed" \
-b "${DIFF}/kidney_Het.bed" \
-u > "${DIFF}/DMR_kidney_heterochromatin.bed"

bedtools intersect \
-a "${DIFF}/mystery_diffmeth_sorted.bed" \
-b "${DIFF}/liver_Het.bed" \
-u > "${DIFF}/DMR_liver_heterochromatin.bed"


echo "Kidney DMRs overlapping heterochromatin:"
wc -l "${DIFF}/DMR_kidney_heterochromatin.bed"

echo "Liver DMRs overlapping heterochromatin:"
wc -l "${DIFF}/DMR_liver_heterochromatin.bed"

# ------------------------------------------------------------
# 4. DMR overlap with differential ATAC peaks
# ------------------------------------------------------------

KIDNEY_ATAC="/vol/COMPEPIWS/groups/shared/ATAC-seq/atacseq3/differential/kidney_accessible_peaks.bed"

LIVER_ATAC="/vol/COMPEPIWS/groups/shared/ATAC-seq/atacseq3/differential/liver_accessible_peaks.bed"


bedtools intersect \
-a "${DIFF}/mystery_diffmeth_sorted.bed" \
-b "${KIDNEY_ATAC}" \
-u > "${DIFF}/DMR_kidney_ATAC.bed"

bedtools intersect \
-a "${DIFF}/mystery_diffmeth_sorted.bed" \
-b "${LIVER_ATAC}" \
-u > "${DIFF}/DMR_liver_ATAC.bed"


echo "Kidney DMRs overlapping accessible ATAC regions:"
wc -l "${DIFF}/DMR_kidney_ATAC.bed"

echo "Liver DMRs overlapping accessible ATAC regions:"
wc -l "${DIFF}/DMR_liver_ATAC.bed"

# ------------------------------------------------------------
# 5. Prepare TSS-centered regions
# ------------------------------------------------------------

awk -vOFS='\t' '
{
    if ($6 == "+") {
        s = $2 - 500
        e = $2 + 1500
    }
    else {
        s = $3 - 1500
        e = $3 + 500
    }

    if (s < 0)
        s = 0

    print $1,s,e
}
' mm10_reduced_chr18_chr19_genes.bed |
sort -k1,1 -k2,2n \
> promoters.bed


# ------------------------------------------------------------
# 6. TSS-centered signal summarization
#
# The resulting matrices can be visualized with:
#   deepTools plotHeatmap
#   deepTools plotProfile
#
# Signals:
#   DNA methylation
#   H3K4me3
#   RNA expression
# ------------------------------------------------------------

# Example structure:
#
# computeMatrix reference-point \
#     --referencePoint TSS \
#     -b 500 \
#     -a 1500 \
#     -R promoters.bed \
#     -S methylation.bedGraph \
#        H3K4me3.bw \
#        RNA.bw \
#     --skipZeros \
#     -o tss_matrix.gz
#
# plotHeatmap \
#     -m tss_matrix.gz \
#     -out tss_heatmap.png
#
# plotProfile \
#     -m tss_matrix.gz \
#     -out tss_profile.png


echo "Integrative analysis preparation completed."
