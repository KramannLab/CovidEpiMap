# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de


#---- Pseudotime analysis on T cell sub-populations with Slingshot and tradeSeq

#devtools::install_github("statOmics/tradeSeq", ref = "fastFitting")
#devtools::install_github("kstreet13/bioc2020trajectories")
library(tradeSeq)
library(Seurat)
library(slingshot)
library(ggplot2)
library(viridis)
library(ggbeeswarm)
library(ggthemes)
library(patchwork)
library(SingleCellExperiment)
set.seed(42)
'%ni%' = Negate('%in%')
source('sc_source/sc_source.R')
indir = '~/sciebo/CovidEpiMap/integrated/'
outdir = '~/sciebo/CovidEpiMap/trajectory_analysis/'


sc = readRDS(file = paste0(indir, 'integrated.RNA.Tcells.annotated.rds'))
DefaultAssay(sc) = 'RNA'
Idents(sc) = 'integrated_annotations'



# Exclude the below cell types as these have a separate 
# developmental origin from the remaining T cell populations
excluded.cell.types = c('Gamma Delta T cells', 'MAIT cells', 
						'CD8+ CD73+ regulatory T cells', 'Atypical NKT cells')

cell.types = names(cell.type.colors)
cell.types = cell.types[cell.types %ni% excluded.cell.types]

sc.subset = subset(sc, integrated_annotations %in% cell.types)
Idents(sc.subset) = factor(Idents(sc.subset), levels = cell.types)



#---- Run slingshot on UMAP embedding

# Use CD8+ naive T cells as progenitor state
start.clus = 'CD8+ naive T cells'
reduction = 'umap'
sds = slingshot(Embeddings(sc.subset, reduction), clusterLabels = Idents(sc.subset), start.clus = start.clus)
sc.subset@tools[['slingshot']] = SlingshotDataSet(sds)



# Plot slingshot curves
pseudotime = slingPseudotime(sds)
curves = colnames(pseudotime)
palette = viridis(100, end = 0.95)


pdf(file = paste0(outdir, 'integrated_Tcells_slingshot_curves.pdf'), width = 12)
plot(reducedDim(sds), col = cell.type.colors[as.vector(sc.subset$integrated_annotations)], pch = 16, cex = 0.5)
lines(sds, lwd = 2, col = 'black')
# Show which curve cells belong to
par(mfrow = c(2, 3))
for (curve in curves) {
  colors = palette[cut(pseudotime[,curve], breaks = 100)]
  print(plot(reducedDim(sds), col = colors, pch = 16, cex = 0.5, main = curve) +
  lines(sds, lwd = 2, col = 'black'))
}
dev.off()



# Add pseudotimes (arclength) to meta data for visualisation
pseudotime2 = slingPseudotime(sds, na = FALSE)
sc.subset$slingshot_pseudotime = pseudotime2[,2]
saveRDS(sc.subset, file = paste0(indir, 'integrated.RNA.Tcells.pseudotime.subset.rds'))


pdf(file = paste0(outdir, 'integrated_Tcells_slingshot_pseudotime.pdf'))
FeaturePlot(sc.subset, feature = 'slingshot_pseudotime', cols = c('red', 'yellow')) + NoLegend()
dev.off()



#---- Plot cells ordered on pseudotime
cell.table = data.frame(cell = colnames(sc.subset), pseudotime = sc.subset$slingshot_pseudotime,
						cell.type = sc.subset$integrated_annotations)


pdf(file = paste0(outdir, 'integrated_Tcells_slingshot_pseudotime_violin.pdf'))
ggplot(cell.table, aes(x = pseudotime, y = cell.type, colour = cell.type)) +
    geom_quasirandom(groupOnX = FALSE) +
    theme_classic() +
    scale_color_manual(values = cell.type.colors, name = 'Cell type') +
    xlab('Slingshot pseudotime (arclength)') +
    ylab('')
dev.off()



#---- Plot condition density along pseudotime

sc.subset$slingshot_pseudotime_curve1 = pseudotime[,1]
sc.subset$slingshot_pseudotime_curve2 = pseudotime[,2]
df = sc.subset@meta.data


pdf(file = paste0(outdir, 'integrated_Tcells_condition_lineage_density.pdf'))
# Lineage 1
ggplot() + 
geom_density(data = df[df$condition %ni% 'healthy',], 
            aes(slingshot_pseudotime_curve1, group = condition, 
                fill = condition), 
            alpha = 0.5, 
            adjust = 2) +
scale_fill_manual(values = viridis(4), name = 'Condition') +
theme_classic() +
xlab('Pseudotime (Lineage 1)') +
ylab('Density')

# Lineage 2
ggplot() + 
geom_density(data = df[df$condition %ni% 'healthy',], 
            aes(slingshot_pseudotime_curve2, group = condition, 
                fill = condition), 
            alpha = 0.5, 
            adjust = 2) +
scale_fill_manual(values = viridis(4), name = 'Condition') +
theme_classic() +
xlab('Pseudotime (Lineage 2)') +
ylab('Density')
dev.off()



pdf(file = paste0(outdir, 'integrated_Tcells_condition_collapsed_lineage_density.pdf'))
# Lineage 1
ggplot() + 
geom_density(data = df[df$condition %ni% 'healthy',], 
            aes(slingshot_pseudotime_curve1, group = condition_collapsed, 
                fill = condition_collapsed), 
            alpha = 0.5, 
            adjust = 2) +
scale_fill_manual(values = viridis(2), name = 'Condition') +
theme_classic() +
xlab('Pseudotime (Lineage 1)') +
ylab('Density')

# Lineage 2
ggplot() + 
geom_density(data = df[df$condition %ni% 'healthy',], 
            aes(slingshot_pseudotime_curve2, group = condition_collapsed, 
                fill = condition_collapsed), 
            alpha = 0.5, 
            adjust = 2) +
scale_fill_manual(values = viridis(2), name = 'Condition') +
theme_classic() +
xlab('Pseudotime (Lineage 2)') +
ylab('Density')
dev.off()



#--- Check for imbalance in local distribution of cells according to condition

scores = bioc2020trajectories::imbalance_score(
  rd = reducedDims(sds),
  cl = sc.subset$condition,
  k = 20, smooth = 40)

grad = viridis::plasma(10, begin = 0, end = 1)
names(grad) = levels(cut(scores$scaled_scores, breaks = 10))

pdf(file = paste0(outdir, 'integrated_Tcells_condition_cell_distribution.pdf'))
plot(reducedDims(sds), 
        col = grad[cut(scores$scaled_scores, breaks = 10)],
        asp = 1, pch = 16, cex = .8, 
        xlab = 'UMAP-1', ylab = 'UMAP-2')
legend('topright', legend = names(grad), 
        col = grad, pch = 16, bty = 'n', cex = 2 / 3)
dev.off()



#---- Fit NB-GAM on gene expression data using tradeSeq

# Fit negative binomial generalized additive model (NB-GAM) on active infection subset

subset = subset(sc.subset, condition %in% c('active_mild', 'active_severe'))
subset = FindVariableFeatures(subset, nfeatures = 10000, verbose = FALSE)
var.genes = VariableFeatures(subset)
subset = subset(subset, features = var.genes)



# As tradeSeq relies on a negative binomial count distribution, use raw counts
counts = as.matrix(subset[['RNA']]@counts)
pseudotime = slingPseudotime(sds, na = FALSE)[colnames(counts),]
cell.weights = slingCurveWeights(sds)[colnames(counts),]



# Fit smoothed average gene expression profile along pseudotime using a NB-GAM
# with a condition-specific smoother for each lineage (run on cluster, took 2.5 days)
BPPARAM = BiocParallel::bpparam()
BPPARAM$workers = 4

sce = fitGAM(counts = counts, pseudotime = pseudotime,
            cellWeights = cell.weights, nknots = 6,
            conditions = as.factor(unname(subset$condition)),
            parallel = TRUE, 
            BPPARAM = BPPARAM)

saveRDS(sce, file = 'sce.parallel.active.subset.rds')


