# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de


#---- Morisita horn overlap of TCR repertoire

library(Seurat)
library(dplyr)
library(stringr)
library(cowplot)
library(vegan)
library(reshape2)
'%ni%' = Negate('%in%')
indir = '~/sciebo/CovidEpiMap/integrated/'
outdir = '~/sciebo/CovidEpiMap/tcr/'
source('sc_source/sc_source.R')


# Read data
sc = readRDS(file = paste0(indir, 'integrated.RNA.Tcells.annotated.rds'))
DefaultAssay(sc) = 'RNA'
Idents(sc) = 'integrated_annotations'

# Subset to MHC class I/II restricted T cells
sc = subset(sc, integrated_annotations %ni% c('Gamma Delta T cells', 'MAIT cells', 'Atypical NKT cells'))

# Compute morisita overlap index of clonotypes between cell types
# The clonotype call is based on the nucleotide sequences. As the chance of overlap between
# NT TCR sequences between patients is very small, we just use he clonotype id 
dm = 1 - as.matrix(vegdist(table(sc$integrated_annotations, sc$Clonotype), method = 'horn'))
coef_matrix = melt(get_upper_triangle(dm, diag = FALSE), na.rm = TRUE)

# Plot
pdf(file = paste0(outdir, 'morisita_overlap_index.pdf'), height = 5, width = 6)
upper_trig_tile_plot(coef_matrix) +
  labs(fill = 'Morisita Horn similarity index')
dev.off()



#---- Plot clonotype expansion groups on UMAP

sc$clonotype_cut = cut(sc$clonotype_size, breaks = c(0, 1, 5, 20, 100, Inf), 
                       labels = rev(c('Hyperexpanded (100 < x <= Inf)', 
                                      'Large (20 < x <= 100)', 
                                      'Medium (5 < x <= 20)', 
                                      'Small (1 < x <= 5)', 
                                      'Single (0 < x <= 1)')))

# Order of groups for plot
group.order = c(NA, rev(c('Hyperexpanded (100 < x <= Inf)', 
                          'Large (20 < x <= 100)', 
                          'Medium (5 < x <= 20)', 
                          'Small (1 < x <= 5)', 
                          'Single (0 < x <= 1)')))

pdf(file = paste0(outdir, 'clonal_expansion_groups.pdf'), width = 9.8)
DimPlot(sc, group.by = 'clonotype_cut', 
        order = TRUE,
        cols = viridis(5))
dev.off()



#---- Bar charts of clonotype expansion groups

df = sc@meta.data
df[is.na(df$Clonotype),'clonotype_cut'] = NA


# Plot
pdf(file = paste0(outdir, 'clonotype_expansion_group_abundance_cell_type.pdf'), width = 6, height = 4)
ggplot(df, aes(x = integrated_annotations, 
               fill = factor(clonotype_cut, levels = group.order, exclude = NULL))) + 
  geom_bar(position = 'fill') + 
  theme_cowplot() + 
  theme(axis.text.x = element_text(angle = 90, size = 8, vjust = 1, hjust = 1, color = 'black'),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = 'black'),
        axis.ticks = element_blank()) +
  scale_fill_viridis(discrete = TRUE, na.value = 'lightgrey') +
  labs(fill = 'Clonal expansion',
       y =  'Relative abundance')
dev.off()

pdf(file = paste0(outdir, 'clonotype_expansion_group_abundance_condition.pdf'), width = 5, height = 4)
ggplot(df, aes(x = condition,
               fill = factor(clonotype_cut, levels = group.order, exclude = NULL))) + 
  geom_bar(position = 'fill') + 
  theme_cowplot() + 
  theme(axis.text.x = element_text(angle = 90, size = 8, vjust = 1, hjust = 1, color = 'black'),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = 'black'),
        axis.ticks = element_blank()) +
  scale_fill_viridis(discrete = TRUE, na.value = 'lightgrey') +
  labs(fill = 'Clonal expansion',
       y =  'Relative abundance')
dev.off()


df$integrated_annotations_condition = factor(df$integrated_annotations_condition,
                                             levels = paste0(rep(names(cell.type.colors),each = 5), 
                                                             levels(df$condition)))

pdf(file = paste0(outdir, 'clonotype_expansion_group_abundance_cell_type_condition.pdf'), width = 10, height = 6)
ggplot(df, aes(x = integrated_annotations_condition, 
               fill = factor(clonotype_cut, levels = group.order, exclude = NULL))) + 
  geom_bar(position = 'fill') + 
  theme_cowplot() + 
  theme(axis.text.x = element_text(angle = 90, size = 8, vjust = 1, hjust = 1, color = 'black'),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = 'black'),
        axis.ticks = element_blank()) +
  scale_fill_viridis(discrete = TRUE, na.value = 'lightgrey') +
  labs(fill = 'Clonal expansion',
       y =  'Relative abundance')
dev.off()
