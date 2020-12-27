# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de


#---- Analysis of dextramer-epitope binding cells 

library(Seurat)
library(dplyr)
library(msigdbr)
library(fgsea)
indir = '~/sciebo/CovidEpiMap/integrated/'
outdir = '~/sciebo/CovidEpiMap/epitope_analysis/'
source('sc_source/sc_source.R')



# GSEA function
run_gsea = function(object, markers, category, plot.title = NULL, subcategory = NULL, out.dir = '.', file.prefix, n = 30){
	# Fetch geneset
	geneSets = msigdbr(species = 'Homo sapiens', category = category, subcategory = subcategory)
	geneSets = geneSets[geneSets$human_gene_symbol %in% rownames(object),]
	m_list = geneSets %>% split(x = .$human_gene_symbol, f = .$gs_name)

	# Run GSEA
	stats = markers$avg_log2FC
	names(stats) = rownames(markers)
	gsea = fgsea(pathways = m_list, stats = stats, minSize = 10, eps = 0.0)
	order = order(gsea$padj, decreasing = FALSE)

	# Plot
	file.name = paste0(out.dir, '/', file.prefix, '_gsea_', paste0(c(category, subcategory), collapse = '_'))

	pdf(file = paste0(file.name, '.pdf'), width = 10, height = 9)
	print(plot_go(gsea.res = gsea,
			gsea.res.order = order, 
			n = n, 
			plot.title = plot.title))
	dev.off()

	# Write to file
	write.table(gsea[order, -8], 
			file = paste0(file.name, '.txt'), 
			sep = '\t', 
			row.names = FALSE, 
			quote = FALSE)
}



#---- Plot dextramer binding on UMAP

sc = readRDS(file = paste0(indir, 'integrated.RNA.Tcells.annotated.rds'))
DefaultAssay(sc) = 'RNA'
Idents(sc) = 'integrated_annotations'


dextramers = c('A0201_5', 'A1101_30', 'A0201_6',
				'A0101_2', 'A1101_29', 'A1101_23',
				'A0201_4', 'A0201_12')

pdf(file = paste0(outdir, 'integrated_Tcells_dextramer_binding.pdf'))
for (dextramer in dextramers){
	print(DimPlot(sc, group.by = dextramer, cols = rev(viridis(2)), order = TRUE) +
		ggtitle(dextramer))
}
dev.off()


# Plot unique A0101-2 binding cells
dex.subset = subset(sc, A0201_5 == 'NO' & A1101_30 == 'NO' & A0201_6 == 'NO' & 
					A0101_2 == 'YES' & A1101_29 == 'NO' & A1101_23 == 'NO' & 
					A0201_4 == 'NO' & A0201_12 == 'NO')

pdf(file = paste0(outdir, 'integrated_Tcells_dextramer_A0101-2_unique_binding.pdf'))
DimPlot(sc, cells.highlight = colnames(dex.subset),
		cols.highlight = viridis(1), sizes.highlight = 0.1)
dev.off()



#---- Plot clonal expansion of non-binding cells

subset = subset(sc, A0201_5 == 'NO' & A1101_30 == 'NO' & A0201_6 == 'NO' & 
					A0101_2 == 'NO' & A1101_29 == 'NO' & A1101_23 == 'NO' & 
					A0201_4 == 'NO' & A0201_12 == 'NO')

pdf(file = paste0(outdir, 'clonal_expansion_non_binding_cells.pdf'))
FeaturePlot(subset, feature = 'clonotype_size', cols = c('yellow', 'red'), order = TRUE)
dev.off()



#---- Plot clonal expansion in A0101-2 binding cells

pdf(file = paste0(outdir, 'clonal_expansion_A0101-2_binding_cells.pdf'))
FeaturePlot(dex.subset, feature = 'clonotype_size', cols = c('yellow', 'red'), order = TRUE)
dev.off()


pdf(file = paste0(outdir, 'clonal_expansion_A0101-2_binding_cells_covid_split.pdf'), width = 10)
FeaturePlot(dex.subset, feature = 'clonotype_size', cols = c('yellow', 'red'), 
			split.by = 'COVID', order = TRUE)
dev.off()



#---- Analysis of A0101-2 binding COVID vs NON-COVID cells

Idents(dex.subset) = 'COVID'

markers = FindMarkers(dex.subset, ident.1 = 'COVID',
      				ident.2 = 'NON-COVID', 
      				min.pct = 0.25,
      				logfc.threshold = 0)

markers = markers[order(markers$avg_log2FC, decreasing = TRUE),]


# Write to file
markers.out = markers[markers$p_val_adj < 0.05,]
markers.out$gene = rownames(markers.out)

write.table(markers.out[,c(6,1:5)], 
        file = paste0(outdir, 'COVID_vs_NON-COVID_A0101-2_binding_cells_dge.txt'),
        quote = FALSE, 
        sep = '\t', 
        row.names = FALSE)


# GSEA
# GO
file.prefix = paste0('COVID_vs_NON-COVID_A0101-2_binding_cells')

run_gsea(object = dex.subset, markers = markers, category = 'C5', subcategory = 'BP',
		out.dir = outdir, plot.title = 'GO',
		file.prefix = file.prefix, n = 40)

# PID
run_gsea(object = dex.subset, markers = markers, category = 'C2', subcategory = 'PID',
		out.dir = outdir, plot.title = 'PID',
		file.prefix = file.prefix)

# Immunological signature
run_gsea(object = dex.subset, markers = markers, category = 'C7',
		out.dir = outdir, plot.title = 'Immunological Signature',
		file.prefix = file.prefix)



#---- A0101-2 binding vs non-binding cells (active severe)

# Only evaluate cell types with >200 binding cells in each group and
# representation from at least 2 patients
dextramer = 'A0101_2'
unique_binders = colnames(dex.subset)
df = subset@meta.data

cell.types = c('CD8+ TEMRA cells', 'CD8+ effector memory T cells 1',
				'CD8+ effector memory T cells 2')


for (cell.type in cell.types){
	# DGEA
	subset = subset(sc, condition == 'active_severe' & integrated_annotations == cell.type)
	Idents(subset) = dextramer

	markers = FindMarkers(subset, ident.1 = 'YES', ident.2 = 'NO', logfc.threshold = 0)
	markers = markers[order(markers$avg_log2FC, decreasing = TRUE),]


	# Write to file
	markers.out = markers[markers$p_val_adj < 0.05,]
	markers.out$gene = rownames(markers.out)
	cell.type.name = gsub(' ', '_', cell.type)

	write.table(markers.out[,c(6, 1:5)], 
			file = paste0(outdir, dextramer, '/active_severe_', cell.type.name, '_binding_vs_nonbinding_dge.txt'), 
			sep = '\t',
			row.names = FALSE, 
			quote = FALSE)


	# GSEA
	# GO
	run_gsea(object = sc, markers = markers, category = 'C5', subcategory = 'BP',
			out.dir = paste0(outdir, dextramer), cell.type = cell.type, 
			file.prefix = paste0('active_severe_', cell.type.name, '_binding_vs_nonbinding'))

	# PID
	run_gsea(object = sc, markers = markers, category = 'C2', subcategory = 'PID',
			out.dir = paste0(outdir, dextramer), cell.type = cell.type, 
			file.prefix = paste0('active_severe_', cell.type.name, '_binding_vs_nonbinding'))

	# Immunological signature
	run_gsea(object = sc, markers = markers, category = 'C7',
			out.dir = paste0(outdir, dextramer), cell.type = cell.type,
			file.prefix = paste0('active_severe_', cell.type.name, '_binding_vs_nonbinding'))	
}



#---- A0101_2-binding recovered severe vs recovered mild CD8+ TEMRA cells 

# Recovered mild binding CD8+ TEMRA cells: 599
# Recovered severe binding CD8+ TEMRA cells: 250


# DGEA
cell.type = 'CD8+ TEMRA cells'
subset = subset(sc, condition %in% c('recovered_mild', 'recovered_severe') & integrated_annotations == cell.type & !!sym(dextramer) == 'YES')
Idents(subset) = 'condition'

markers = FindMarkers(subset, ident.1 = 'recovered_severe', ident.2 = 'recovered_mild', logfc.threshold = 0)
markers = markers[order(markers$avg_log2FC, decreasing = TRUE),]


# Write to file
markers.out = markers[markers$p_val_adj < 0.05,]
markers.out$gene = rownames(markers.out)
cell.type.name = gsub(' ', '_', cell.type)
file.name = paste0(cell.type.name, '_binding_recovered_severe_vs_recovered_mild')

write.table(markers.out[,c(6, 1:5)], 
			file = paste0(outdir, dextramer, '/', file.name, '_dge.txt'), 
			sep = '\t',
			row.names = FALSE, 
			quote = FALSE)


# GSEA
# GO
run_gsea(object = sc, markers = markers, category = 'C5', subcategory = 'BP',
		out.dir = paste0(outdir, dextramer), cell.type = cell.type, 
		file.prefix = file.name)

# PID
run_gsea(object = sc, markers = markers, category = 'C2', subcategory = 'PID',
		out.dir = paste0(outdir, dextramer), cell.type = cell.type, 
		file.prefix = file.name)

# Immunological signature
run_gsea(object = sc, markers = markers, category = 'C7',
		out.dir = paste0(outdir, dextramer), cell.type = cell.type,
		file.prefix = file.name)	









