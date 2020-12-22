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
library(msigdbr)
library(gridExtra)
library(fgsea)
library(slingshot)
library(knitr)
library(UpSetR)
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



#---- Test for DE genes along pseudotime

association.res = associationTest(sce, lineages = TRUE, l2fc = log2(2))
association.res$padj_lineage1_conditionactive_mild = p.adjust(association.res$pvalue_lineage1_conditionactive_mild, 'fdr')
association.res$padj_lineage2_conditionactive_mild = p.adjust(association.res$pvalue_lineage2_conditionactive_mild, 'fdr')
association.res$padj_lineage1_conditionactive_severe = p.adjust(association.res$pvalue_lineage1_conditionactive_severe, 'fdr')
association.res$padj_lineage2_conditionactive_severe = p.adjust(association.res$pvalue_lineage2_conditionactive_severe, 'fdr')
association.res$gene = rownames(association.res)


# Write to file
write.table(association.res[,c(21,1:20)], 
			file = paste0(outdir, 'association.test.active.subset.txt'),
			sep = '\t', 
			row.names = FALSE,
			quote = FALSE)


# Order and subset data
order.mild.1 = order(association.res$waldStat_lineage1_conditionactive_mild, decreasing = TRUE)
order.mild.2 = order(association.res$waldStat_lineage2_conditionactive_mild, decreasing = TRUE)
order.severe.1 = order(association.res$waldStat_lineage1_conditionactive_severe, decreasing = TRUE)
order.severe.2 = order(association.res$waldStat_lineage2_conditionactive_severe, decreasing = TRUE)

assoc.order.mild.1 = association.res[order.mild.1,]
assoc.order.mild.2 = association.res[order.mild.2,]
assoc.order.severe.1 = association.res[order.severe.1,]
assoc.order.severe.2 = association.res[order.severe.2,]

mild.genes.1 = rownames(assoc.order.mild.1)[which(assoc.order.mild.1$padj_lineage1_conditionactive_mild < 0.05)]
mild.genes.2 = rownames(assoc.order.mild.2)[which(assoc.order.mild.2$padj_lineage2_conditionactive_mild < 0.05)]
severe.genes.1 = rownames(assoc.order.severe.1)[which(assoc.order.severe.1$padj_lineage1_conditionactive_severe < 0.05)]
severe.genes.2 = rownames(assoc.order.severe.2)[which(assoc.order.severe.2$padj_lineage2_conditionactive_severe < 0.05)]


# Plot top ranking genes
plot_top_genes(gene.set = mild.genes.1, model = sce, n = 30, 
				out.dir = outdir, file.name = 'association_test_lineage1_top_genes_mild')

plot_top_genes(gene.set = mild.genes.2, model = sce, n = 30, 
				out.dir = outdir, file.name = 'association_test_lineage2_top_genes_mild')

plot_top_genes(gene.set = severe.genes.1, model = sce, n = 30, 
				out.dir = outdir, file.name = 'association_test_lineage1_top_genes_severe')

plot_top_genes(gene.set = severe.genes.2, model = sce, n = 30, 
				out.dir = outdir, file.name = 'association_test_lineage2_top_genes_severe')


# Plot overlap of significant genes
pdf(file = paste0(outdir, 'association_test_intersection_genes.pdf'))
upset(fromList(list('Active mild (lineage 1)' = mild.genes.1, 'Active mild (lineage 2)' = mild.genes.2)))
upset(fromList(list('Active severe (lineage 1)' = severe.genes.1, 'Active severe (lineage 2)' = severe.genes.2)))
upset(fromList(list('Active mild (lineage 1)' = mild.genes.1, 'Active severe (lineage 1)' = severe.genes.1)))
upset(fromList(list('Active mild (lineage 2)' = mild.genes.2, 'Active severe (lineage 2)' = severe.genes.2)))
dev.off()



#---- Test for DE between start- and end-point of lineages (progenitor markers)

start.res = startVsEndTest(sce, l2fc = log2(2))
start.res$padj = p.adjust(start.res$pvalue, 'fdr')

order = order(start.res$waldStat, decreasing = TRUE)
start.order = start.res[order,]
start.order$gene = rownames(start.order)
start.genes = rownames(start.order)[which(start.order$padj < 0.05)]


# Write to file
write.table(start.order[,c(7,1:6)], 
			file = paste0(outdir, 'start.test.active.subset.txt'),
			sep = '\t', 
			row.names = FALSE,
			quote = FALSE)


# Plot top ranking genes
plot_top_genes(gene.set = start.genes, model = sce, n = 30, 
				out.dir = outdir, file.name = 'start_test_top_genes')



#---- Test for DE between end-points of lineages (differentiated markers)

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
plot_top_genes(gene.set = end.genes, model = sce, n = length(end.genes), 
				out.dir = outdir, file.name = 'diff_end_test_top_genes')



#---- Test DE where the lineages split

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
plot_top_genes(gene.set = split.genes, model = sce, n = length(split.genes), 
				out.dir = outdir, file.name = 'lineage_split_test_top_genes')



#---- Test for DE genes between conditions

condition.res = conditionTest(sce, l2fc = log2(2), global = TRUE, pairwise = TRUE)
condition.res$padj_lineage1 = p.adjust(condition.res$pvalue_lineage1, 'fdr')
condition.res$padj_lineage2 = p.adjust(condition.res$pvalue_lineage2, 'fdr')
condition.res$gene = rownames(condition.res)


# Write to file
write.table(condition.res[,c(12,1:11)], 
			file = paste0(outdir, 'condition.test.active.subset.txt'),
			sep = '\t', 
			row.names = FALSE,
			quote = FALSE)


# Order and subset data
order.1 = order(condition.res$waldStat_lineage1, decreasing = TRUE)
order.2 = order(condition.res$waldStat_lineage2, decreasing = TRUE)
cond.order.1 = condition.res[order.1,]
cond.order.2 = condition.res[order.2,]
genes.1 = rownames(cond.order.1)[which(cond.order.1$padj_lineage1 < 0.05)]
genes.2 = rownames(cond.order.2)[which(cond.order.2$padj_lineage1 < 0.05)]


# Plot top ranking genes
plot_top_genes(gene.set = genes.1, model = sce, n = 100, 
				out.dir = outdir, file.name = 'condition_test_lineage1_top_genes')

plot_top_genes(gene.set = genes.2, model = sce, n = 100, 
				out.dir = outdir, file.name = 'condition_test_lineage2_top_genes')


# Heatmap of DE genes between conditions
# Lineage 1
smooth = predictSmooth(sce, gene = genes.1, nPoints = 100, tidy = FALSE)

pdf(file = paste0(outdir, 'condition_test_lineage1_DEG_heatmap.pdf'), width = 6)
plot_heatmaps(smooth.res = smooth, 
			cond1.start = 201, cond1.end = 300, 
			cond2.start = 1, cond2.end = 100,
			main1 = 'Active severe (lineage 1)',
			main2 = 'Active mild (lineage 1)')
dev.off()


# Lineage 2
smooth = predictSmooth(sce, gene = genes.2, nPoints = 100, tidy = FALSE)

pdf(file = paste0(outdir, 'condition_test_lineage2_DEG_heatmap.pdf'), width = 6)
plot_heatmaps(smooth.res = smooth, 
			cond1.start = 301, cond1.end = 400, 
			cond2.start = 101, cond2.end = 200,
			main1 = 'Active severe (lineage 2)',
			main2 = 'Active mild (lineage 2)')
dev.off()



#---- GSEA with fgsea

# Fetch genesets
geneSets = msigdbr(species = 'Homo sapiens', category = 'C5', subcategory = 'BP')
geneSets = geneSets[geneSets$human_gene_symbol %in% names(sce),]
m_list = geneSets %>% split(x = .$human_gene_symbol, f = .$gs_name)



#--- GSEA of condition test results

stats.1 = cond.order.1$waldStat_lineage1
stats.2 = cond.order.2$waldStat_lineage2
names(stats.1) = rownames(cond.order.1)
names(stats.2) = rownames(cond.order.2)

gsea.1 = fgsea(pathways = m_list, stats = stats.1, minSize = 10, eps = 0.0)
gsea.2 = fgsea(pathways = m_list, stats = stats.2, minSize = 10, eps = 0.0)

gsea.order.1 = order(gsea.1$padj, decreasing = FALSE)
gsea.order.2 = order(gsea.2$padj, decreasing = FALSE)


# Write to file
write.table(gsea.1[gsea.order.1, -8], 
			file = paste0(outdir, 'condition.test.active.subset.gsea.lineage1.txt'), 
			sep = '\t', 
			row.names = FALSE, 
			quote = FALSE)

write.table(gsea.2[gsea.order.2, -8], 
			file = paste0(outdir, 'condition.test.active.subset.gsea.lineage2.txt'), 
			sep = '\t', 
			row.names = FALSE, 
			quote = FALSE)


# Plot GO terms
p1 = plot_go(gsea.res = gsea.1, gsea.res.order = gsea.order.1, n = 10, plot.title = 'Lineage 1')
p1 = plot_go(gsea.res = gsea.2, gsea.res.order = gsea.order.2, n = 10, plot.title = 'Lineage 2')


pdf(file = paste0(outdir, 'condition_test_active_subset_gsea.pdf'), width = 12, height = 4)
cowplot::plot_grid(p1, p2)
dev.off()



#---- GSEA of start- vs end-point test results

stats.start = start.order$waldStat
names(stats.start) = rownames(start.order)
stats.start = stats.start[!is.na(stats.start)]

gsea.start = fgsea(pathways = m_list, stats = stats.start, minSize = 10, eps = 0.0)
gsea.start.order = order(gsea.start$padj, decreasing = FALSE)


# Write to file
write.table(gsea.start[gsea.start.order, -8], 
			file = paste0(outdir, 'start.test.active.subset.gsea.txt'), 
			sep = '\t', 
			row.names = FALSE, 
			quote = FALSE)


# Plot GO terms
pdf(file = paste0(outdir, 'start_test_active_subset_gsea.pdf'), height = 5)
plot_go(gsea.res = gsea.start, gsea.res.order = gsea.start.order)
dev.off()

