# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de


#---- Plot differential expression results

library(dplyr)
source('sc_source/sc_source.R')
indir = '~/sciebo/CovidEpiMap/diff_expression/diff_genes/'

# Genes to plot
genes = read.table(file = paste0(indir, 'analysis_diff_genes/relevant_genes_final_plot_DEG.txt'), header = TRUE)
genes = genes$gene

# Get DGEA results
pattern = '.active_severe.vs.active_mild.txt'
de.files = list.files(path = indir, pattern = pattern)
dge.table = lapply(file.path(paste0(indir,de.files)), read.table, sep = '\t', header = TRUE)

cell.types = sub(pattern, '', de.files)
cell.types = sub('integrated.diff.genes.', '', cell.types)
names(dge.table) = cell.types

pdf(file = paste0(indir, '../active_severe_vs_active_mild_relevant_DEG_final_plot.pdf'), height = 11, width = 5)
plot_dge_nice(dge.table = dge.table, genes = genes)
dev.off()



#---- Epitope A0101-2-binding cells (severe vs mild)

# Get DGEA results
indir = '~/sciebo/CovidEpiMap/epitope_analysis/severe_vs_mild/'
pattern = 'severe_vs_mild_dge.txt'
de.files = list.files(path = indir, pattern = pattern)
dge.table = lapply(file.path(paste0(indir,de.files)), read.table, sep = '\t', header = TRUE)

cell.types = sub(pattern, '', de.files)
cell.types = sub('_A0101-2_binding_', '', cell.types)
names(dge.table) = cell.types

# Genes to plot
genes = read.table(file = paste0(indir, 'DEG_epitope_specific_cells.txt'), header = TRUE)
genes = genes$gene

pdf(file = paste0(indir, '../A0101-2_binding_severe_vs_mild_selected_diff_genes.pdf'), height = 12.2, width = 4)
plot_dge_nice(dge.table = dge.table, genes = genes)
dev.off()

