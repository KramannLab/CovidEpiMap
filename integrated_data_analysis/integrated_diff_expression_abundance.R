# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de


#---- Differential gene expression analysis across conditions of integrated data

library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggrepel)
library(emmeans)
indir = '~/sciebo/CovidEpiMap/integrated/'
source('sc_source/sc_source.R')


sc = readRDS(file = paste0(indir, 'integrated.RNA.Tcells.annotated.rds'))
DefaultAssay(sc) = 'RNA'


cell.types = levels(Idents(sc))
sc$celltype.condition = paste(Idents(sc), sc$condition)
Idents(sc) = 'celltype.condition'



for (cell.type in cell.types){
	# Healthy vs active mild
	condition = 'active_mild'
	control = 'healthy'

	markers = get_markers(sc, condition = condition, 
		control = control, cell.type = cell.type)

	top_heatmap(sc, markers = markers, cell.type = cell.type, 
		disease = condition, control = control, n = 15)


	# Healthy vs active severe
	condition = 'active_severe'
	control = 'healthy'

	markers = get_markers(sc, condition = condition, 
		control = control, cell.type = cell.type)
	
	top_heatmap(sc, markers = markers, cell.type = cell.type, 
		disease = condition, control = control, n = 15)


	# Mild vs severe active
	condition = 'active_severe'
	control = 'active_mild'

	markers = get_markers(sc, condition = condition, 
		control = control, cell.type = cell.type)
	
	top_heatmap(sc, markers = markers, cell.type = cell.type, 
		disease = condition, control = control, n = 15)


	# Mild vs severe recovered
	condition = 'recovered_severe'
	control = 'recovered_mild'

	markers = get_markers(sc, condition = condition, 
		control = control, cell.type = cell.type)
	
	top_heatmap(sc, markers = markers, cell.type = cell.type, 
		disease = condition, control = control, n = 15)
}



#---- Plot relevant diff genes

indir = '~/sciebo/CovidEpiMap/integrated/diff_genes/analysis_diff_genes/'
genes = read.table(file = paste0(indir, 'relevant_genes_plot_DEG.txt'), header = TRUE)
genes = genes$gene.name

pdf(file = paste0(indir, 'relevant_genes_DGE.pdf'), height = 10)
DotPlot(sc, features = rev(genes), group.by = 'condition', dot.scale = 8) + 
coord_flip() + 
scale_colour_gradient2(low = 'blue', mid = 'lightgrey', high = 'red', 
    midpoint = 0, limits = c(-2,2), oob = scales::squish) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off() 


genes = read.table(file = paste0(indir, 'relevant_genes_interactions_plot_DEG.txt'), header = TRUE)
genes = genes$gene


pdf(file = paste0(indir, 'relevant_genes_interactions_DGE.pdf'), height = 12)
for (cell.type in names(cell.type.colors)){
	subset = subset(sc, integrated_annotations == cell.type)

	print(DotPlot(subset, features = rev(genes), group.by = 'condition', dot.scale = 8) + 
	coord_flip() + 
	scale_colour_gradient2(low = 'blue', mid = 'lightgrey', high = 'red', 
					midpoint = 0, limits = c(-2,2), oob = scales::squish) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(cell.type))
}
dev.off()


pdf(file = paste0(indir, 'relevant_genes_interactions_DGE_per_patient.pdf'), height = 12, width = 12)
for (cell.type in names(cell.type.colors)){
	subset = subset(sc, integrated_annotations == cell.type)

	print(DotPlot(subset, features = rev(genes), group.by = 'patient', dot.scale = 8) + 
	coord_flip() + 
	scale_colour_gradient2(low = 'blue', mid = 'lightgrey', high = 'red', 
					midpoint = 0, limits = c(-2,2), oob = scales::squish) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(cell.type))
}
dev.off()



#---- Cell-type abundance testing
# https://www.nxn.se/valent/2020/11/28/s9jjv32ogiplagwx8xrkjk532p7k28
outdir = '~/sciebo/CovidEpiMap/cell_type_abundance/'

# Transform data
df = sc@meta.data


# Count total number of cells per patient and number of cells per cell type
data = df %>%
select(condition, patient, integrated_annotations) %>%
group_by(patient) %>%
count(integrated_annotations, name = 'count', .drop = FALSE) %>%
add_count(wt = count, name = 'total') %>%
mutate(other = total - count) %>%
as.data.frame()


# Add condition to table again
condition = rep(c('healthy', 'active_mild',
				'active_severe', 'recovered_mild',
				'recovered_severe'), 
				each = 3*13)
data = cbind(data, condition)


# Save counts
write.table(data, file = paste0(outdir, 'integrated.cell.type.counts.txt'),
			sep = '\t', row.names = FALSE, quote = FALSE)



# Healthy vs active mild
conditions = c('healthy', 'active_mild')
plots = get_abundance(df = data, condition.vector = conditions)

pdf(file = paste0(outdir, conditions[1], '_vs_', conditions[2], '.pdf'))
plots$p1
plots$p2
dev.off()


# Healthy vs active severe
conditions = c('healthy', 'active_severe')
plots = get_abundance(df = data, condition.vector = conditions)

pdf(file = paste0(outdir, conditions[1], '_vs_', conditions[2], '.pdf'))
plots$p1
plots$p2
dev.off()


# Healthy vs recovered mild
conditions = c('healthy', 'recovered_mild')
plots = get_abundance(df = data, condition.vector = conditions)

pdf(file = paste0(outdir, conditions[1], '_vs_', conditions[2], '.pdf'))
plots$p1
plots$p2
dev.off()


# Healthy vs recovered severe
conditions = c('healthy', 'recovered_severe')
plots = get_abundance(df = data, condition.vector = conditions)

pdf(file = paste0(outdir, conditions[1], '_vs_', conditions[2], '.pdf'))
plots$p1
plots$p2
dev.off()


# Actuve mild vs active severe
conditions = c('active_mild', 'active_severe')
plots = get_abundance(df = data, condition.vector = conditions)

pdf(file = paste0(outdir, conditions[1], '_vs_', conditions[2], '.pdf'))
plots$p1
plots$p2
dev.off()



