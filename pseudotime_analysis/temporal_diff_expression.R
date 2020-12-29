# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de


#---- Temporal differential gene expression analysis with tradeSeq

# Sources
# https://kstreet13.github.io/bioc2020trajectories/articles/workshopTrajectories.html
# https://statomics.github.io/tradeSeq/articles/tradeSeq.html
library(tradeSeq)
library(ggplot2)
library(viridis)
library(Seurat)
library(pheatmap)
library(SingleCellExperiment)
library(gridExtra)
library(slingshot)
set.seed(42)
source('sc_source/sc_source.R')
indir = '~/sciebo/CovidEpiMap/integrated/'
outdir = '~/sciebo/CovidEpiMap/trajectory_analysis/'


# Get slingshot data
sc = readRDS(file = paste0(indir, 'integrated.RNA.Tcells.pseudotime.subset.rds'))
sds = sc@tools$slingshot

# NB-GAM model based on active mild/severe subset
sce = readRDS(file = paste0(outdir, 'sce.parallel.active.subset.rds'))



#---- Visualize knots on UMAP

pdf(file = paste0(outdir, 'tradeseq_knots.pdf'))
plotGeneCount(sds, models = sce, clusters = apply(slingClusterLabels(sds), 1, which.max))
dev.off()



#---- Test for DE between start- and end-point of lineages (progenitor markers)

# Temporal DE
start.res = startVsEndTest(sce, l2fc = log2(2), lineages = TRUE)
start.res$padj_lineage1 = p.adjust(start.res$pvalue_lineage1, 'fdr')
start.res$padj_lineage2 = p.adjust(start.res$pvalue_lineage2, 'fdr')

order = order(start.res$waldStat, decreasing = TRUE)
start.order = start.res[order,]
start.order$gene = rownames(start.order)

# Write to file
write.table(start.order[,c(14,1:13)], 
			file = paste0(outdir, 'start.vs.end.test.active.subset.txt'),
			sep = '\t', 
			row.names = FALSE,
			quote = FALSE)

# Plot top ranking genes
sig.subset.1 = start.res[which(start.res$padj_lineage1 < 0.05),]
sig.subset.2 = start.res[which(start.res$padj_lineage2 < 0.05),]
order.1 = order(sig.subset.1$waldStat_lineage1, decreasing = TRUE)
order.2 = order(sig.subset.2$waldStat_lineage2, decreasing = TRUE)
start.genes.1 = rownames(sig.subset.1[order.1,])
start.genes.2 = rownames(sig.subset.2[order.2,])

file.prefix.1 = 'start_vs_end_test_lineage1'
file.prefix.2 = 'start_vs_end_test_lineage2'

plot_top_genes(gene.set = start.genes.1, model = sce, n = 30, 
				out.dir = outdir, file.name = paste0(file.prefix.1, '_top_genes'))

plot_top_genes(gene.set = start.genes.2, model = sce, n = 30, 
				out.dir = outdir, file.name = paste0(file.prefix.2, '_top_genes'))


# GSEA
bg.genes = names(sce)
order.1 = order(start.res$waldStat_lineage1, decreasing = TRUE)
order.2 = order(start.res$waldStat_lineage2, decreasing = TRUE)

stats.1 = start.res[order.1,]$waldStat_lineage1
stats.2 = start.res[order.2,]$waldStat_lineage2
names(stats.1) = rownames(start.res[order.1,])
names(stats.2) = rownames(start.res[order.2,])
stats.1 = stats.1[!is.na(stats.1)]
stats.2 = stats.2[!is.na(stats.2)]

# GO
run_gsea(bg.genes = bg.genes, stats = stats.1, 
		category = 'C5', subcategory = 'BP',
		out.dir = outdir, plot.title = 'GO',
		file.prefix = file.prefix.1, n = 30)

run_gsea(bg.genes = bg.genes, stats = stats.2, 
		category = 'C5', subcategory = 'BP',
		out.dir = outdir, plot.title = 'GO',
		file.prefix = file.prefix.2, n = 30)

# PID
run_gsea(bg.genes = bg.genes, stats = stats.1, 
		category = 'C2', subcategory = 'PID',
		out.dir = outdir, plot.title = 'PID',
		file.prefix = file.prefix.1, n = 30)

run_gsea(bg.genes = bg.genes, stats = stats.2, 
		category = 'C2', subcategory = 'PID',
		out.dir = outdir, plot.title = 'PID',
		file.prefix = file.prefix.2, n = 30)

# Immunological signature
run_gsea(bg.genes = bg.genes, stats = stats.1, 
		category = 'C7',
		out.dir = outdir, plot.title = 'Immunological Signature',
		file.prefix = file.prefix.1, n = 30)

run_gsea(bg.genes = bg.genes, stats = stats.2, 
		category = 'C7',
		out.dir = outdir, plot.title = 'Immunological Signature',
		file.prefix = file.prefix.2, n = 30)



#---- Test for DE between end-points of lineages (differentiated markers)

# Temporal DE
end.res = diffEndTest(sce, l2fc = log2(2))
end.res$padj = p.adjust(end.res$pvalue, 'fdr')

order = order(end.res$waldStat, decreasing = TRUE)
end.order = end.res[order,]
end.order$gene = rownames(end.order)
end.genes = rownames(end.order)[which(end.order$padj < 0.05)]

# Write to file
write.table(end.order[,c(6,1:5)], 
			file = paste0(outdir, 'diff.end.test.active.subset.txt'),
			sep = '\t', 
			row.names = FALSE,
			quote = FALSE)

# Plot top ranking genes
file.prefix = 'diff_end_test'

plot_top_genes(gene.set = end.genes, model = sce, n = length(end.genes), 
				out.dir = outdir, file.name = paste0(file.prefix, '_top_genes'))


# GSEA
stats = end.order$waldStat
names(stats) = rownames(end.order)
stats = stats[!is.na(stats)]

# GO
run_gsea(bg.genes = bg.genes, stats = stats, 
		category = 'C5', subcategory = 'BP',
		out.dir = outdir, plot.title = 'GO',
		file.prefix = file.prefix, n = 30)

# PID
run_gsea(bg.genes = bg.genes, stats = stats, 
		category = 'C2', subcategory = 'PID',
		out.dir = outdir, plot.title = 'PID',
		file.prefix = file.prefix, n = 30)

# Immunological signature
run_gsea(bg.genes = bg.genes, stats = stats, 
		category = 'C7',
		out.dir = outdir, plot.title = 'Immunological Signature',
		file.prefix = file.prefix, n = 30)



#---- Test DE where the lineages split

# Temporal DE
# Use lower l2fc threshold to capture more subtle changes
split.res = earlyDETest(sce, knots = c(2, 3), l2fc = log2(1.5))
split.res$padj = p.adjust(split.res$pvalue, 'fdr')

order = order(split.res$waldStat, decreasing = TRUE)
split.order = split.res[order,]
split.order$gene = rownames(split.order)
split.genes = rownames(split.order)[which(split.order$padj < 0.05)]

# Write to file
write.table(split.order[,c(6,1:5)], 
			file = paste0(outdir, 'lineage.split.test.active.subset.txt'),
			sep = '\t', 
			row.names = FALSE,
			quote = FALSE)

# Plot top ranking genes
file.prefix = 'lineage_split_test'
plot_top_genes(gene.set = split.genes, model = sce, n = length(split.genes), 
				out.dir = outdir, file.name = paste0(file.prefix, '_top_genes'))


# GSEA
stats = split.order$waldStat
names(stats) = rownames(split.order)
stats = stats[!is.na(stats)]

# GO
run_gsea(bg.genes = bg.genes, stats = stats, 
         category = 'C5', subcategory = 'BP',
         out.dir = outdir, plot.title = 'GO',
         file.prefix = file.prefix, n = 30)

# PID
run_gsea(bg.genes = bg.genes, stats = stats, 
         category = 'C2', subcategory = 'PID',
         out.dir = outdir, plot.title = 'PID',
         file.prefix = file.prefix, n = 30)

# Immunological signature
run_gsea(bg.genes = bg.genes, stats = stats, 
         category = 'C7',
         out.dir = outdir, plot.title = 'Immunological Signature',
         file.prefix = file.prefix, n = 30)



#---- Test for DE genes between conditions

# Temporal DE
condition.res = conditionTest(sce, l2fc = log2(2), global = TRUE, pairwise = TRUE)
condition.res$padj_lineage1 = p.adjust(condition.res$pvalue_lineage1, 'fdr')
condition.res$padj_lineage2 = p.adjust(condition.res$pvalue_lineage2, 'fdr')

order = order(condition.res$waldStat, decreasing = TRUE)
condition.order = condition.res[order,]
condition.order$gene = rownames(condition.order)

# Write to file
write.table(condition.order[,c(12,1:11)], 
            file = paste0(outdir, 'condition.test.active.subset.txt'),
            sep = '\t', 
            row.names = FALSE,
            quote = FALSE)

# Plot top ranking genes
sig.subset.1 = condition.res[which(condition.res$padj_lineage1 < 0.05),]
sig.subset.2 = condition.res[which(condition.res$padj_lineage2 < 0.05),]
order.1 = order(sig.subset.1$waldStat_lineage1, decreasing = TRUE)
order.2 = order(sig.subset.2$waldStat_lineage2, decreasing = TRUE)
genes.1 = rownames(sig.subset.1[order.1,])
genes.2 = rownames(sig.subset.2[order.2,])

file.prefix.1 = 'condition_test_lineage1'
file.prefix.2 = 'condition_test_lineage2'

plot_top_genes(gene.set = genes.1, model = sce, n = 100, 
				out.dir = outdir, file.name = paste0(file.prefix.1, '_top_genes'))

plot_top_genes(gene.set = genes.2, model = sce, n = 100, 
				out.dir = outdir, file.name = paste0(file.prefix.2, '_top_genes'))

# Heatmap
# Lineage 1
smooth = predictSmooth(sce, gene = genes.1, nPoints = 100, tidy = FALSE)

pdf(file = paste0(outdir, file.prefix.1, '_DEG_heatmap.pdf'), width = 6)
plot_heatmaps(smooth.res = smooth, 
			cond1.start = 201, cond1.end = 300, 
			cond2.start = 1, cond2.end = 100,
			main1 = 'Active severe (lineage 1)',
			main2 = 'Active mild (lineage 1)')
dev.off()

# Lineage 2
smooth = predictSmooth(sce, gene = genes.2, nPoints = 100, tidy = FALSE)

pdf(file = paste0(outdir, file.prefix.2, '_DEG_heatmap.pdf'), width = 6)
plot_heatmaps(smooth.res = smooth, 
			cond1.start = 301, cond1.end = 400, 
			cond2.start = 101, cond2.end = 200,
			main1 = 'Active severe (lineage 2)',
			main2 = 'Active mild (lineage 2)')
dev.off()


# GSEA
order.1 = order(condition.res$waldStat_lineage1, decreasing = TRUE)
order.2 = order(condition.res$waldStat_lineage2, decreasing = TRUE)

stats.1 = condition.res[order.1,]$waldStat_lineage1
stats.2 = condition.res[order.2,]$waldStat_lineage2
names(stats.1) = rownames(condition.res[order.1,])
names(stats.2) = rownames(condition.res[order.2,])
stats.1 = stats.1[!is.na(stats.1)]
stats.2 = stats.2[!is.na(stats.2)]

# GO
run_gsea(bg.genes = bg.genes, stats = stats.1, 
         category = 'C5', subcategory = 'BP',
         out.dir = outdir, plot.title = 'GO',
         file.prefix = file.prefix.1, n = 30)

run_gsea(bg.genes = bg.genes, stats = stats.2, 
         category = 'C5', subcategory = 'BP',
         out.dir = outdir, plot.title = 'GO',
         file.prefix = file.prefix.2, n = 30)

# PID
run_gsea(bg.genes = bg.genes, stats = stats.1, 
         category = 'C2', subcategory = 'PID',
         out.dir = outdir, plot.title = 'PID',
         file.prefix = file.prefix.1, n = 30)

run_gsea(bg.genes = bg.genes, stats = stats.2, 
         category = 'C2', subcategory = 'PID',
         out.dir = outdir, plot.title = 'PID',
         file.prefix = file.prefix.2, n = 30)

# Immunological signature
run_gsea(bg.genes = bg.genes, stats = stats.1, 
         category = 'C7',
         out.dir = outdir, plot.title = 'Immunological Signature',
         file.prefix = file.prefix.1, n = 30)

run_gsea(bg.genes = bg.genes, stats = stats.2, 
         category = 'C7',
         out.dir = outdir, plot.title = 'Immunological Signature',
         file.prefix = file.prefix.2, n = 30)

