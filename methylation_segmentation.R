# ============================================================
# Mouse Methylome Analysis
# Methylation-state segmentation
# Input: Tissue-level merged methylomes
# Method: MethylSeekR through RnBeads
# States: UMR, LMR, HMR, PMD
# ============================================================

library(RnBeads)

project_dir <- "/vol/COMPEPIWS/groups/wgbs3/tasks/results"

report_dir <- file.path(
    project_dir,
    "rnbeads_report"
)

segmentation_dir <- file.path(
    project_dir,
    "segmentation"
)

segmentation_bed_dir <- file.path(
    project_dir,
    "segmentation_bed"
)

dir.create(
    segmentation_dir,
    showWarnings = FALSE
)

dir.create(
    segmentation_bed_dir,
    showWarnings = FALSE
)

# ------------------------------------------------------------
# Load merged methylome
# ------------------------------------------------------------

rnb.merged <- load.rnb.set(
    file.path(
        report_dir,
        "rnbSet_with_mystery"
    )
)

# ------------------------------------------------------------
# Segment kidney methylome
# ------------------------------------------------------------

rnb.merged <- rnb.execute.segmentation(
    rnb.set = rnb.merged,
    sample.name = "kidney",
    chr.sel = "chr18",
    n.cores = 1,
    plot.path = segmentation_dir
)

# ------------------------------------------------------------
# Segment liver methylome
# ------------------------------------------------------------

rnb.merged <- rnb.execute.segmentation(
    rnb.set = rnb.merged,
    sample.name = "liver",
    chr.sel = "chr18",
    n.cores = 1,
    plot.path = segmentation_dir
)

# ------------------------------------------------------------
# Export kidney segmentation
# ------------------------------------------------------------

rnb.bed.from.segmentation(
    rnb.set = rnb.merged,
    sample.name = "kidney",
    type = "final",
    store.path = segmentation_bed_dir
)

# ------------------------------------------------------------
# Export liver segmentation
# ------------------------------------------------------------

rnb.bed.from.segmentation(
    rnb.set = rnb.merged,
    sample.name = "liver",
    type = "final",
    store.path = segmentation_bed_dir
)

# ------------------------------------------------------------
# Save segmented object
# ------------------------------------------------------------

save.rnb.set(
    rnb.merged,
    file.path(
        report_dir,
        "rnbSet_segmented"
    ),
    archive = FALSE
)

cat(
    "Methylation segmentation completed.\n"
)
