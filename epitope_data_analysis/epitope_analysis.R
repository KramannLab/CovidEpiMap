# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de


#---- Analysis of dextramer-epitope binding cells 

library(Seurat)
library(dplyr)
library(ggplot2)
library(viridis)
library(viper)
library(tidyverse)
library(progeny)
library(dorothea)
'%ni%' = Negate('%in%')
indir = '~/sciebo/CovidEpiMap/integrated/'
outdir = '~/sciebo/CovidEpiMap/epitope_analysis/'
source('sc_source/sc_source.R')


sc = readRDS(file = paste0(indir, 'integrated.RNA.Tcells.annotated.rds'))
DefaultAssay(sc) = 'RNA'
Idents(sc) = 'integrated_annotations'

# Subset to MHC class I/II restricted T cells
sc = subset(sc, integrated_annotations %ni% c('Gamma Delta T cells', 'MAIT cells', 'Atypical NKT cells'))



#---- Plot dextramer binding on UMAP

dextramers = c('A0201_5', 'A1101_30', 'A0201_6',
				'A0101_2', 'A1101_29', 'A1101_23',
				'A0201_4', 'A0201_12')

pdf(file = paste0(outdir, 'integrated_Tcells_dextramer_binding.pdf'))
for (dextramer in dextramers){
	print(DimPlot(sc, group.by = dextramer, cols = rev(viridis(2)), order = TRUE) +
		ggtitle(dextramer))
}
dev.off()



#---- Plot unique A0101-2 binding cells

dex.subset = subset(sc, A0201_5 == 'NO' & A1101_30 == 'NO' & A0201_6 == 'NO' & 
					A0101_2 == 'YES' & A1101_29 == 'NO' & A1101_23 == 'NO' & 
					A0201_4 == 'NO' & A0201_12 == 'NO')

saveRDS(dex.subset, file = paste0(indir, 'A0101-2.unique.binders.rds'))

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



#---- Clonotype size of A0101-2 binding cells vs pseudotime (severe and mild)

# Get pseudotimes
sc.subset = readRDS(file = paste0(indir, 'integrated.RNA.Tcells.pseudotime.subset.rds'))
df = sc.subset@meta.data
df = df[df$condition %ni% 'healthy',]
df$integrated_annotations = factor(df$integrated_annotations, 
                                    levels = names(cell.type.colors))

# Subset to A0101-2 binding cells
unique.binders = colnames(dex.subset)
df = df[rownames(df) %in% unique.binders,]


# Scatter plot
pdf(file = paste0(outdir, 'A0101-2_binding_cells_clonotype_size_vs_pseudotime.pdf'), width = 7, height = 4)
# Lineage 1
ggplot(df, aes(x = slingshot_pseudotime_curve1, y = clonotype_size)) + 
geom_point(aes(colour = integrated_annotations)) + 
scale_colour_manual(values = cell.type.colors) +
stat_smooth(method = 'loess', colour = 'black') +
facet_wrap(. ~ condition_collapsed) +
theme_classic() +
theme(strip.background = element_rect(colour = 'black', fill = 'lightgrey')) +
xlab('Pseudotime (Lineage 1)') +
ylab('Clonotype size')
# Lineage 2
ggplot(df, aes(x = slingshot_pseudotime_curve2, y = clonotype_size)) + 
geom_point(aes(colour = integrated_annotations)) + 
scale_colour_manual(values = cell.type.colors) +
stat_smooth(method = 'loess', colour = 'black') +
facet_wrap(. ~ condition_collapsed) +
theme_classic() +
theme(strip.background = element_rect(colour = 'black', fill = 'lightgrey')) +
xlab('Pseudotime (Lineage 2)') +
ylab('Clonotype size')
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

file.prefix = 'COVID_vs_NON-COVID_A0101-2_binding_cells'

write.table(markers.out[,c(6,1:5)], 
        file = paste0(outdir, file.prefix, '_dge.txt'),
        quote = FALSE, 
        sep = '\t', 
        row.names = FALSE)


# GSEA
bg.genes = rownames(dex.subset)
stats = markers$avg_log2FC
names(stats) = rownames(markers)
stats = stats[!is.na(stats)]

# GO
run_gsea(bg.genes = bg.genes, stats = stats, 
		category = 'C5', subcategory = 'BP',
		out.dir = outdir, plot.title = 'GO',
		file.prefix = file.prefix, n = 40)

# PID
run_gsea(bg.genes = bg.genes, stats = stats, 
		category = 'C2', subcategory = 'PID',
		out.dir = outdir, plot.title = 'PID',
		file.prefix = file.prefix)

# Immunological signature
run_gsea(bg.genes = bg.genes, stats = stats, 
		category = 'C7',
		out.dir = outdir, plot.title = 'Immunological Signature',
		file.prefix = file.prefix)



#---- A0101-2 binding vs non-binding cells (active severe)

# Only evaluate cell types with >200 binding cells in each group and
# representation from at least 2 patients
cell.types = c('CD8+ TEMRA cells', 'CD8+ effector memory T cells 1',
				'CD8+ effector memory T cells 2')

for (cell.type in cell.types){
	binding = colnames(subset(dex.subset, condition == 'active_severe' & integrated_annotations == cell.type))
	non.binding = colnames(subset(sc, condition == 'active_severe' & integrated_annotations == cell.type & A0101_2 == 'NO'))
	subset = subset(sc, cells = c(binding, non.binding))
	Idents(subset) = 'A0101_2'


	# DGEA
	markers = FindMarkers(subset, ident.1 = 'YES', ident.2 = 'NO', logfc.threshold = 0)
	markers = markers[order(markers$avg_log2FC, decreasing = TRUE),]

	# Write to file
	markers.out = markers[markers$p_val_adj < 0.05,]
	markers.out$gene = rownames(markers.out)
	cell.type.name = gsub(' ', '_', cell.type)

	file.prefix = paste0('active_severe_', cell.type.name, '_binding_vs_nonbinding')
	write.table(markers.out[,c(6, 1:5)], 
			file = paste0(outdir, file.prefix, '_dge.txt'), 
			sep = '\t',
			row.names = FALSE, 
			quote = FALSE)


	# GSEA
	stats = markers$avg_log2FC
	names(stats) = rownames(markers)
	stats = stats[!is.na(stats)]

	# GO
	run_gsea(bg.genes = bg.genes, stats = stats, 
			category = 'C5', subcategory = 'BP',
			out.dir = outdir, plot.title = 'GO',
			file.prefix = file.prefix, n = 40)

	# PID
	run_gsea(bg.genes = bg.genes, stats = stats, 
			category = 'C2', subcategory = 'PID',
			out.dir = outdir, plot.title = 'PID',
			file.prefix = file.prefix)

	# Immunological signature
	run_gsea(bg.genes = bg.genes, stats = stats, 
			category = 'C7',
			out.dir = outdir, plot.title = 'Immunological Signature',
			file.prefix = file.prefix)
}



#---- A0101-2 binding recovered severe vs recovered mild

cell.types = c('CD8+ TEMRA cells', 'CD8+ effector memory T cells 1')

for (cell.type in cell.types){
	subset = subset(dex.subset, integrated_annotations == cell.type)
	Idents(subset) = 'condition'

	# DGEA
	markers = FindMarkers(subset, ident.1 = 'recovered_severe', 
						ident.2 = 'recovered_mild', logfc.threshold = 0)
	markers = markers[order(markers$avg_log2FC, decreasing = TRUE),]

	# Write to file
	markers.out = markers[markers$p_val_adj < 0.05,]
	markers.out$gene = rownames(markers.out)
	cell.type.name = gsub(' ', '_', cell.type)

	file.prefix = paste0(cell.type.name, '_A0101-2_binding_recovered_severe_vs_recovered_mild')
	write.table(markers.out[,c(6, 1:5)], 
			file = paste0(outdir, file.prefix, '_dge.txt'), 
			sep = '\t',
			row.names = FALSE, 
			quote = FALSE)


	# GSEA
	stats = markers$avg_log2FC
	names(stats) = rownames(markers)
	stats = stats[!is.na(stats)]

	# GO
	run_gsea(bg.genes = bg.genes, stats = stats, 
			category = 'C5', subcategory = 'BP',
			out.dir = outdir, plot.title = 'GO',
			file.prefix = file.prefix, n = 40)

	# PID
	run_gsea(bg.genes = bg.genes, stats = stats, 
			category = 'C2', subcategory = 'PID',
			out.dir = outdir, plot.title = 'PID',
			file.prefix = file.prefix)

	# Immunological signature
	run_gsea(bg.genes = bg.genes, stats = stats, 
			category = 'C7',
			out.dir = outdir, plot.title = 'Immunological Signature',
			file.prefix = file.prefix)
}



#---- A0101-2 binding severe vs mild

cell.types = c('CD8+ TEMRA cells', 'CD8+ effector memory T cells 1')

for (cell.type in cell.types){
	subset = subset(dex.subset, integrated_annotations == cell.type)
	Idents(subset) = 'condition_collapsed'

	# DGEA
	markers = FindMarkers(subset, ident.1 = 'severe', 
						ident.2 = 'mild', logfc.threshold = 0)
	markers = markers[order(markers$avg_log2FC, decreasing = TRUE),]

	# Write to file
	markers.out = markers[markers$p_val_adj < 0.05,]
	markers.out$gene = rownames(markers.out)
	cell.type.name = gsub(' ', '_', cell.type)
	
	file.prefix = paste0(cell.type.name, '_A0101-2_binding_severe_vs_mild')
	write.table(markers.out[,c(6, 1:5)], 
			file = paste0(outdir, file.prefix, '_dge.txt'), 
			sep = '\t',
			row.names = FALSE, 
			quote = FALSE)


	# GSEA
	stats = markers$avg_log2FC
	names(stats) = rownames(markers)
	stats = stats[!is.na(stats)]

	# GO
	run_gsea(bg.genes = bg.genes, stats = stats, 
			category = 'C5', subcategory = 'BP',
			out.dir = outdir, plot.title = 'GO',
			file.prefix = file.prefix, n = 40)

	# PID
	run_gsea(bg.genes = bg.genes, stats = stats, 
			category = 'C2', subcategory = 'PID',
			out.dir = outdir, plot.title = 'PID',
			file.prefix = file.prefix)

	# Immunological signature
	run_gsea(bg.genes = bg.genes, stats = stats, 
			category = 'C7',
			out.dir = outdir, plot.title = 'Immunological Signature',
			file.prefix = file.prefix)
}

