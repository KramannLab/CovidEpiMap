# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de


#---- Differential gene expression analysis across conditions of integrated data

library(Seurat)
library(dplyr)
library(ggplot2)
indir = '~/sciebo/CovidEpiMap/integrated/'
source('sc_source/sc_source.R')


sc = readRDS(file = paste0(indir, 'integrated.RNA.Tcells.annotated.rds'))
DefaultAssay(sc) = 'RNA'


cell.types = levels(Idents(sc))
sc$celltype.condition = paste(Idents(sc), sc$condition)
Idents(sc) = 'celltype.condition'


# Cell type vs cell type across conditions
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


# Global DE between COVID vs NON-COVID
outdir = '~/sciebo/CovidEpiMap/diff_expression/'
Idents(sc) = 'COVID'
condition = 'COVID'
control = 'NON-COVID'

markers = FindMarkers(sc, ident.1 = condition, 
                      ident.2 = control, 
                      min.pct = 0.25, 
                      logfc.threshold = 0)
markers$gene = rownames(markers)

write.table(markers[,c(6, 1:5)], 
           file = paste0(outdir, 'integrated.diff.genes.', condition, '.vs.', control, '.txt'),
           quote = FALSE, 
           sep = '\t', 
           row.names = FALSE)



#---- Plot relevant diff genes

indir = '~/sciebo/CovidEpiMap/diff_expression/diff_genes/analysis_diff_genes/'
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



#---- Plot differential expression of selected genes

# Source
# https://davemcg.github.io/post/lets-plot-scrna-dotplots/

# Genes to plot
genes = read.table(file = paste0(indir, 'relevant_genes_final_plot_DEG.txt'), header = TRUE)
genes = genes$gene

# Get DGEA results
indir = '~/sciebo/CovidEpiMap/diff_expression/diff_genes/'
pattern = '.active_severe.vs.active_mild.txt'
de.files = list.files(path = indir, pattern = pattern)
dge.table = lapply(file.path(paste0(indir,de.files)), read.table, sep = '\t', header = TRUE)

cell.types = sub(pattern, '', de.files)
cell.types = sub('integrated.diff.genes.', '', cell.types)
names(dge.table) = cell.types

gene.list = list()
for (i in 1:length(dge.table)){
  table = dge.table[[i]]
  
  col = table %>% 
        filter(gene %in% genes) %>% 
        select(gene, avg_log2FC, p_val_adj) %>% 
        mutate(cluster = names(dge.table[i]),
               log10_padj = -log10(p_val_adj))
  
  gene.list[[i]] = col
}
gene.table = do.call(rbind, gene.list)
gene.table$log10_padj_edit = gene.table$log10_padj
gene.table[gene.table$log10_padj > 15,'log10_padj_edit'] = 15

pdf(file = paste0(indir, 'active_severe_vs_active_mild_relevant_DEG_final_plot.pdf'), height = 11, width = 5)
gene.table %>%
ggplot(aes(x = cluster, y = gene, color = avg_log2FC, size = log10_padj_edit)) +
geom_point() +
scale_color_viridis_c(name = 'Log2FC') + 
cowplot::theme_cowplot() + 
theme(axis.line  = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      axis.ticks = element_blank()) +
labs(x = '',
     y = '',
     size = '-log10(pAdj)')
dev.off()

