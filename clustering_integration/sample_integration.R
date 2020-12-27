# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de


#---- Integration of all samples with Harmony

library(Seurat)
library(rlist)
library(dplyr)
library(ggplot2)
library(harmony)
options(future.globals.maxSize = 30720*1024^2)
outdir = '~/sciebo/CovidEpiMap/integrated/'
source('sc_source/sc_source.R')
'%ni%' = Negate('%in%')


# Prepare meta data table
patients = c('1', '2', '3', 
			'5', '6', '7', 
			'11', '12', '13', 
			'15', '18', '19', 
			'29', '31', '32')
condition = c('active_mild', 'active_mild', 'active_severe', 
			'recovered_mild', 'recovered_mild', 'recovered_mild', 
			'active_severe', 'recovered_severe', 'recovered_severe',
			'active_mild', 'recovered_severe', 'active_severe',
			'healthy', 'healthy', 'healthy')
condition.table = as.data.frame(cbind(patients, condition))
rownames(condition.table) = paste0(patients, '_GEX_SURF')



# Integration based on GEX data only

assay = 'RNA'
obj.list = c()
for (patient in patients){
	sample = paste0(patient, '_GEX_SURF')
	indir = paste0('~/sciebo/CovidEpiMap/per_patient/', sample, '/')
	sc = readRDS(file = paste0(indir, sample, '.rds'))

	# Find variable features again
	DefaultAssay(sc) = assay
	sc = FindVariableFeatures(sc, nfeatures = 2000, verbose = FALSE)
	sc = RenameCells(object = sc, add.cell.id = sample)
	obj.list = list.append(obj.list, sc)
}


# Integrate and re-cluster
sc = merge(obj.list[[1]], obj.list[-1])
sc = NormalizeData(sc, verbose = FALSE)
sc = FindVariableFeatures(sc, verbose = FALSE)
sc = ScaleData(sc, features = rownames(sc), verbose = FALSE)
sc = RunPCA(sc, verbose = FALSE)
sc = RunHarmony(sc, group.by.vars = 'orig.ident')
sc = RunUMAP(sc, reduction = 'harmony', dims = 1:30, verbose = FALSE)
sc = FindNeighbors(sc, reduction = 'harmony', dims = 1:30, verbose = FALSE)
sc = FindClusters(sc, resolution = 1, verbose = FALSE)


# Add meta data 
sc$patient = sub('_GEX_SURF', '', sc$orig.ident)
sc$condition = condition.table[sc$orig.ident, 'condition']



#---- Visualize counts from ADT library (dextramer + CITE seq)
DefaultAssay(sc) = 'ADT'
sc = ScaleData(sc, verbose = FALSE, features = rownames(sc))
markers = rownames(sc[['ADT']])


pdf(file = paste0(outdir, 'integrated_RNA_ADT_Markers_Heatmap.pdf'), width = 15)
DoHeatmap(subset(sc, downsample = 2000), features = markers, assay = 'ADT', 
	raster = FALSE, size = 1) + NoLegend()
dev.off()



#---- Remove poor quality clusters and non-T cells

bad.clusters = c('5', '6', '8', '11', '12', 
				'16', '20', '21', '22', '23', '24')
sc.subset = subset(sc, subset = RNA_snn_res.1 %ni% bad.clusters)



#---- Re-integrate and re-cluster the T cells only

sc.subset = NormalizeData(sc.subset, verbose = FALSE)
sc.subset = FindVariableFeatures(sc.subset, verbose = FALSE)
sc.subset = ScaleData(sc.subset, features = rownames(sc.subset), verbose = FALSE)
sc.subset = RunPCA(sc.subset, verbose = FALSE)
sc.subset = RunHarmony(sc.subset, group.by.vars = 'orig.ident')
sc.subset = RunUMAP(sc.subset, reduction = 'harmony', dims = 1:30, verbose = FALSE)
sc.subset = FindNeighbors(sc.subset, reduction = 'harmony', dims = 1:30, verbose = FALSE)
sc.subset = FindClusters(sc.subset, resolution = 0.5, verbose = FALSE)


# Remove poor quality cluster 
bad.clusters = c('7')
sc.subset = subset(sc.subset, subset = RNA_snn_res.0.5 %ni% bad.clusters)

sc.subset = RenameIdents(sc.subset, `0` = 'CD8+ TEMRA cells', `1` = 'CD8+ effector memory T cells 1',
					`2` = 'CD8+ central memory T cells', `3` = 'CD8+ naive T cells',
					`4`  = 'CD8+ effector memory T cells 2', `5` = 'CD8+ NK-like TEMRA cells',
					`6` = 'MAIT cells', `8` = 'CD8+ NK-like early effector T cells',
					`9` = 'CD8+ CD73+ regulatory T cells', `10` = 'Gamma Delta T cells',
					`11` = 'CD8+ cycling effector T cells', `12` = 'Atypical NKT cells',
					`13` = 'CD8+ exhausted T cells')


# Order meta data
sc.subset$integrated_annotations = Idents(sc.subset)

sc.subset$integrated_annotations = factor(sc.subset$integrated_annotations, 
									levels = names(cell.type.colors))

sc.subset$condition = factor(sc.subset$condition,
						levels = c('healthy', 'active_mild', 
									'active_severe', 'recovered_mild',
									'recovered_severe'))

sc.subset$patient = factor(sc.subset$patient,
						levels = c('29', '31', '32',
									'1', '2', '15',
									'3', '11', '19',
									'5', '6', '7',
									'12', '13', '18'))



#---- Visualize UMAP
pdf(file = 'integrated_RNA_Tcells.pdf', width = 10)
DimPlot(sc.subset, cols = cell.type.colors)
dev.off()



#---- Bar chart with TCR V genes in each cluster (show top 20 per cluster)

Idents(sc.subset) = sc.subset$integrated_annotations


pdf(file = paste0(outdir, 'plots/integrated_Tcells_TCRV_genes.pdf'), width = 12)
for (cell.type in names(cell.type.colors)){
	subset = subset(sc.subset, idents = cell.type)

	# Get all TCR V genes
	cell.table = data.frame(cell = colnames(subset), cluster = Idents(subset),
						tcrv = subset$TCR_V_GENE)

	# Get top 20 TCR V genes
	plot.table = cell.table %>% 
	group_by(tcrv) %>% 
	count() %>% 
	arrange(desc(n)) %>% 
	as.data.frame %>% head(20)

	print(ggplot(cell.table %>% filter(tcrv %in% plot.table$tcrv), 
				aes(x = tcrv, fill = as.factor(cluster))) + 
		geom_bar() + 
		labs(y = 'Count', x = element_blank(), fill = 'Cell type') +
		ggtitle(paste0(cell.type, ' (n = ', nrow(cell.table), ')')) +
		theme_classic() +
		theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
		scale_fill_manual(values = cell.type.colors[cell.type]))
}
dev.off()



#---- Correlate cell annotations 

sc.subset$integrated_annotations_condition = paste0(sc.subset$integrated_annotations, sc.subset$condition)
Idents(sc.subset) = 'integrated_annotations_condition'

# RNA
condition.vector = levels(sc.subset$condition)
p = correlation_heatmap(object = sc.subset, assay = 'RNA', conditionVector = condition.vector)

pdf(file = paste0(outdir, 'celltype_condition_correlation.pdf'), width = 15, height = 10)
draw(p)
dev.off()


# ADT
p = correlation_heatmap(object = sc.subset, assay = 'ADT', conditionVector = condition.vector)

pdf(file = paste0(outdir, 'celltype_condition_correlation_ADT.pdf'), width = 15, height = 10)
draw(p)
dev.off()



#--- Cell-cycle analysis

sc.subset = CellCycleScoring(sc.subset, s.features = cc.genes$s.genes, 
							g2m.features = cc.genes$g2m.genes)

pdf(file = paste0(outdir, 'integrated_Tcells_cell_cycle_score.pdf'), width = 12)
FeaturePlot(sc.subset, features = c('S.Score', 'G2M.Score'))
VlnPlot(sc.subset, features = 'S.Score', pt.size = 0, split.by = 'condition', sort = TRUE)
VlnPlot(sc.subset, features = 'G2M.Score', pt.size = 0, split.by = 'condition', sort = TRUE)
VlnPlot(sc.subset, features = 'S.Score', pt.size = 0, split.by = 'condition_condensed', sort = TRUE)
VlnPlot(sc.subset, features = 'G2M.Score', pt.size = 0, split.by = 'condition_condensed', sort = TRUE)
dev.off()



#---- Percent mt mapping across conditions

sc.subset[['percent.ribo']] = PercentageFeatureSet(sc.subset, pattern = '^RP[SL]')
sc.subset[['percent.mt']] = PercentageFeatureSet(sc.subset, pattern = '^MT-')


pdf(file = paste0(outdir, 'integrated_Tcells_percent_mt_condition.pdf'))
VlnPlot(sc.subset, feature = 'percent.mt', pt.size = 0, group.by = 'condition')
VlnPlot(sc.subset, feature = 'percent.mt', pt.size = 0, group.by = 'condition_condensed')
dev.off()



#---- Save data

saveRDS(sc.subset, file = 'integrated.RNA.Tcells.annotated.rds')


