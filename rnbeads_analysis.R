# Mouse Methylome Analysis — RnBeads
# Developing mouse kidney and liver, E14.5/E15.5; 8 single-end WGBS samples; mm10/GRCm38 reduced reference.

library(RnBeads); library(ggplot2); library(matrixStats); library(reshape2); library(GenomicRanges)

# Paths
project_dir <- "/vol/COMPEPIWS/groups/wgbs3/tasks/results"
methylation_dir <- file.path(project_dir, "bismark/methylation_calls/methylation_coverage")
report_dir <- file.path(project_dir, "rnbeads_report")
dir.create(report_dir, showWarnings=FALSE, recursive=TRUE)

# Sample annotation
sample_annotation <- file.path(project_dir, "sample_annotation.csv")
annotation <- data.frame(
  sampleId=c("kidney_14.5.1","kidney_14.5.2","kidney_15.5.1","kidney_15.5.2","liver_14.5.1","liver_14.5.2","liver_15.5.1","liver_15.5.2"),
  tissue=c(rep("kidney",4),rep("liver",4)),
  timepoint=c("14.5","14.5","15.5","15.5","14.5","14.5","15.5","15.5"),
  filename=c("kidney_14.5.1_trimmed_bismark_bt2.deduplicated.bismark.cov.gz","kidney_14.5.2_trimmed_bismark_bt2.deduplicated.bismark.cov.gz","kidney_15.5.1_trimmed_bismark_bt2.deduplicated.bismark.cov.gz","kidney_15.5.2_trimmed_bismark_bt2.deduplicated.bismark.cov.gz","liver_14.5.1_trimmed_bismark_bt2.deduplicated.bismark.cov.gz","liver_14.5.2_trimmed_bismark_bt2.deduplicated.bismark.cov.gz","liver_15.5.1_trimmed_bismark_bt2.deduplicated.bismark.cov.gz","liver_15.5.2_trimmed_bismark_bt2.deduplicated.bismark.cov.gz")
)
write.csv(annotation, sample_annotation, row.names=FALSE)

# RnBeads configuration
rnb.options(assembly="mm10", identifiers.column="sampleId", differential.comparison.columns=c("tissue","timepoint"), differential.report.sites=FALSE, import.bed.style="bismarkCov", filtering.coverage.threshold=10, filtering.missing.value.quantile=1, filtering.high.coverage.outliers=TRUE, filtering.low.coverage.masking=TRUE, import.table.separator=",")

# Import methylation data
data_source <- c(methylation_dir, sample_annotation, "filename")
rnb.set <- rnb.execute.import(data.source=data_source, data.type="bs.bed.dir")
print(rnb.set)

# Run complete RnBeads analysis
result <- rnb.run.analysis(dir.reports=report_dir, data.source=data_source, data.type="bs.bed.dir")

# Load preprocessed methylation data
rnb.set <- load.rnb.set(file.path(report_dir,"rnbSet_preprocessed"))
print(rnb.set)

# PCA
meth.mat <- meth(rnb.set)
vars <- apply(meth.mat,1,var,na.rm=TRUE)
meth.mat.complete <- meth.mat[!is.na(vars) & vars>0,]
meth.mat.complete <- meth.mat.complete[complete.cases(meth.mat.complete),]
pca <- prcomp(t(meth.mat.complete),scale.=TRUE)

pca_scores <- as.data.frame(pca$x[,1:2]); pca_scores$Sample <- rownames(pca_scores)
pca_scores$Tissue <- pheno(rnb.set)$tissue[match(pca_scores$Sample,pheno(rnb.set)$sampleId)]
pca_scores$Timepoint <- pheno(rnb.set)$timepoint[match(pca_scores$Sample,pheno(rnb.set)$sampleId)]

pca_plot <- ggplot(pca_scores,aes(PC1,PC2,label=Sample,colour=Tissue)) + geom_point(size=4) + geom_text(vjust=-0.8,size=3) + theme_classic() + labs(title="PCA of DNA methylation profiles",x="PC1",y="PC2")
ggsave(file.path(project_dir,"PCA_methylation.png"),pca_plot,width=8,height=6)

# PCA variance explained
variance_explained <- pca$sdev^2/sum(pca$sdev^2)
scree_df <- data.frame(PC=seq_along(variance_explained),Variance=variance_explained)
scree_plot <- ggplot(scree_df,aes(PC,Variance)) + geom_point() + geom_line() + theme_classic() + labs(title="Variance explained by principal components",x="Principal component",y="Proportion of variance")
ggsave(file.path(project_dir,"PCA_screeplot.png"),scree_plot,width=7,height=5)

# Methylation density by tissue
meth_long <- as.data.frame(meth.mat); meth_long$CpG <- rownames(meth_long)
meth_long <- melt(meth_long,id.vars="CpG",variable.name="Sample",value.name="Methylation")
meth_long <- merge(meth_long,pheno(rnb.set),by.x="Sample",by.y="sampleId")

density_plot <- ggplot(meth_long,aes(Methylation,colour=tissue)) + geom_density(na.rm=TRUE) + theme_classic() + labs(title="Methylation distribution by tissue",x="DNA methylation",y="Density")
ggsave(file.path(project_dir,"methylation_density_tissue.png"),density_plot,width=8,height=6)

# CpG-context methylation
context_plot <- rnb.plot.density(rnb.set,region.type="cpgislands")
ggsave(file.path(project_dir,"methylation_density_cpg_context.png"),context_plot,width=8,height=6)

# Top 1,000 variable CpGs
meth.complete <- meth.mat[complete.cases(meth.mat),]
row.sd <- rowSds(meth.complete)
top1000 <- order(row.sd,decreasing=TRUE)[1:1000]
heat.mat <- meth.complete[top1000,]
write.csv(heat.mat,file.path(project_dir,"top1000_variable_CpGs.csv"))

# Custom mystery-region annotation
mystery_annotation <- "/vol/COMPEPIWS/data/annotation/annotation_mm10_mystery.RData"
rnb.load.annotation(mystery_annotation,"mystery")
rnb.set <- summarize.regions(rnb.set,"mystery")
save.rnb.set(rnb.set,file.path(report_dir,"rnbSet_with_mystery"),archive=FALSE)

# Top 100 variable mystery regions
mystery.mat <- meth(rnb.set,"mystery")
mystery.complete <- mystery.mat[complete.cases(mystery.mat),]
mystery.sd <- rowSds(mystery.complete)
top100_mystery <- order(mystery.sd,decreasing=TRUE)[1:100]
heat.mystery <- mystery.complete[top100_mystery,]
write.csv(heat.mystery,file.path(project_dir,"top100_variable_mystery_regions.csv"))

# Merge samples by tissue
rnb.merged <- mergeSamples(rnb.set,"tissue")
print(rnb.merged)

# Export tissue-level methylation tracks
bedgraph_dir <- file.path(project_dir,"bedGraph")
dir.create(bedgraph_dir,showWarnings=FALSE)
rnb.RnBSet.to.bedGraph(rnb.merged,bedgraph_dir)

# Differential methylation
dm <- rnb.execute.diffMeth(rnb.merged,comparison="tissue")
print(dm)

# Save analysis objects
saveRDS(rnb.set,file.path(project_dir,"rnbSet_preprocessed.rds"))
saveRDS(rnb.merged,file.path(project_dir,"rnbSet_tissue_merged.rds"))

cat("RnBeads methylome analysis completed.\n")
