# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de


#---- Violin plots of canonical cell type markers

library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(genesorteR)
library(viridis)
source('sc_source/sc_source.R')


indir = '~/sciebo/CovidEpiMap/integrated/'
sc = readRDS(file = paste0(indir, 'integrated.RNA.Tcells.annotated.rds'))
DefaultAssay(sc) = 'RNA'
Idents(sc) = 'integrated_annotations'

# Genes to plot
genes = read.table(file = paste0(indir, 'marker_expression/selected_cell_type_markers.txt'), header = TRUE)
genes = genes$gene

p = get_violin(object = sc, features.use = genes)
pdf(file = paste0(indir, 'cell_type_markers.pdf'), width = 8, height = 4)
p
dev.off()

get_violin(object = sc, features.use = genes[1:3])

#---- Bar charts with cell type distribution

cell.table = data.frame(cell = colnames(sc), condition = sc$condition,
                        cluster = sc$integrated_annotations, 
                        patient = sc$patient)


pdf(file = paste0(outdir, 'integrated_Tcells_barchart_celltype_condition.pdf'), width = 8)
ggplot(cell.table, aes(x = condition, fill = as.factor(cluster))) + 
geom_bar(position = position_fill(reverse = TRUE)) + 
labs(y = 'Proportion', x = element_blank(), fill = 'Cell type') +
theme_classic() +
theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
scale_fill_manual(values = cell.type.colors)
dev.off()

pdf(file = paste0(outdir, 'integrated_Tcells_barchart_celltype_patient.pdf'), width = 20)
ggplot(cell.table, aes(x = patient, fill = as.factor(cluster))) + 
geom_bar(position = position_fill(reverse = TRUE)) + 
labs(y = 'Proportion', x = element_blank(), fill = 'Cell type') +
theme_classic() +
theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
scale_fill_manual(values = cell.type.colors)
dev.off()



#---- Bar charts with average cell type distribution

df = cell.table %>%
group_by(patient) %>%
count(cluster, name = 'count', .drop = FALSE) %>%
add_count(wt = count, name = 'total') %>%
mutate(fraction = count/total) %>%
as.data.frame

condition = rep(c('healthy', 'active_mild',
        'active_severe', 'recovered_mild',
        'recovered_severe'), 
        each = 3*13)
df = cbind(df, condition)

df = df %>% 
group_by(condition, cluster) %>%
summarize(mean = mean(fraction)) %>%
as.data.frame

df$condition = factor(df$condition, levels = c('healthy', 'active_mild',
        'active_severe', 'recovered_mild',
        'recovered_severe'))
df$condition_collapsed = sub('active_', '', df$condition)
df$condition_collapsed = sub('recovered_', '', df$condition_collapsed)


pdf(file = paste0(indir, 'integrated_Tcells_barchart_celltype_condition_average.pdf'), width = 8)
ggplot(df, aes(fill = cluster, y = mean, x = condition)) + 
    geom_bar(stat = 'identity', position = position_fill(reverse = TRUE)) + 
labs(y = 'Average proportion', x = element_blank(), fill = 'Cell type') +
theme_classic() +
theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
scale_fill_manual(values = cell.type.colors)
dev.off()

pdf(file = paste0(indir, 'integrated_Tcells_barchart_celltype_condition_collapsed_average.pdf'))
ggplot(df, aes(fill = cluster, y = mean, x = condition_collapsed)) + 
    geom_bar(stat = 'identity', position = position_fill(reverse = TRUE)) + 
labs(y = 'Average proportion', x = element_blank(), fill = 'Cell type') +
theme_classic() +
theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
scale_fill_manual(values = cell.type.colors)
dev.off()



#---- DotPlot of top10 genes based on specificity score per cluster

sg = run_genesorter(sc, assay = 'RNA', slot = 'data', write.file = TRUE,
  out.dir = indir, file.name.prefix = 'integrated.tcells.')
specScore = sg$specScore

genes = list()
for (i in 1:ncol(specScore)){
  genes[[i]] = rownames(head(specScore[order(specScore[,i], decreasing = TRUE),], 10))
}

top.genes = unique(unlist(genes))
height = length(top.genes) * 0.2
width = length(levels(Idents(sc))) * 0.5

pdf(file = paste0(indir, 'integrated_Tcells_GeneSorterTop10.pdf'), height = height, width = width)
DotPlot(sc, features = rev(top.genes), assay = 'RNA') + 
coord_flip() + 
theme(axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.x = element_blank(),
      axis.title.y = element_blank())
dev.off()



#---- Plot expression of immune genes

genes = c('IL2', 'IL6', 'IL10', 'IL12A', 
          'IL12B', 'IL13', 'IL15', 'TNF')

pdf(file = paste0(indir, 'integrated_Tcells_immune_markers.pdf'), width = 8)
DotPlot(sc, features = rev(genes), assay = 'RNA', group.by = 'condition') + 
coord_flip() + 
theme(axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.x = element_blank(),
      axis.title.y = element_blank())

DotPlot(sc, features = rev(genes), assay = 'RNA', group.by = 'integrated_annotations') + 
coord_flip() + 
theme(axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.x = element_blank(),
      axis.title.y = element_blank())
dev.off()

pdf(file = paste0(indir, 'integrated_Tcells_immune_markers_violin.pdf'), width = 15)
for (gene in genes){
  print(VlnPlot(sc, feature = gene, 
                split.by = 'condition', 
                group.by = 'integrated_annotations'))
}
dev.off()



#---- Dotplot of activation markers

genes = c('TNFRSF9', 'CD69', 'CD44',
            'SELL', 'HLA-DR', 'CD38',
            'CD40LG', 'KLRG1', 'LAMP1',
            'XCL1', 'XCL2', 'CX3CR1')


pdf(file = paste0(indir, 'integrated_Tcells_activation_markers.pdf'))
for (cell.type in names(cell.type.colors)){
  subset = subset(sc, integrated_annotations == cell.type)

  print(DotPlot(subset, features = rev(genes), group.by = 'condition', dot.scale = 8) + 
  coord_flip() + 
  scale_colour_gradient2(low = 'blue', mid = 'lightgrey', high = 'red', 
          midpoint = 0, limits = c(-2,2), oob = scales::squish) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    ggtitle(cell.type))
}
dev.off()


# Activation marker expression in A0101-2 binding cells

outdir = '~/sciebo/CovidEpiMap/epitope_analysis/'
unique.binders = colnames(subset(sc, A0101_2 == 'YES' & A0201_5 == 'NO' & A1101_30 == 'NO' &
         A1101_23 == 'NO' & A1101_29 == 'NO' & A0201_4 == 'NO' & 
         A0201_6 == 'NO' & A0201_12 == 'NO'))

sc[['A0101_2_binding']] = ifelse(colnames(sc) %in% unique.binders, 'unique_binder', 'non_binder')
cell.types = c('CD8+ TEMRA cells', 'CD8+ effector memory T cells 1', 'CD8+ effector memory T cells 2')

pdf(file = paste0(outdir, 'active_severe_A0101_2_unique_binding_non_binding_activation_markers.pdf'))
for (cell.type in cell.types){
  subset = subset(sc, integrated_annotations == cell.type & condition == 'active_severe')
  subset[['ADT_UPDATED']] = NULL
  
  print(DotPlot(subset, features = rev(genes), group.by = 'A0101_2_binding', dot.scale = 8, assay = 'RNA') + 
          coord_flip() + 
          scale_colour_gradient2(low = 'blue', mid = 'lightgrey', high = 'red', 
                                 midpoint = 0, limits = c(-2,2), oob = scales::squish) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1),
                axis.title.x = element_blank(),
                axis.title.y = element_blank()) +
          ggtitle(cell.type))
}
dev.off()


# Activation marker subset expression in A0101-2 binding cells
genes = c('CD69', 'CD44', 'HLA-DR', 'CD38')
cell.types = c('CD8+ TEMRA cells', 'CD8+ effector memory T cells 1')

sc$condition_collapsed_A0101_2_binding = paste(sc$condition_collapsed, sc$A0101_2_binding, sep = '_')
Idents(sc) = 'condition_collapsed_A0101_2_binding'

pdf(file = paste0(outdir, 'A0101_2_unique_binding_non_binding_activation_markers.pdf'))
for (cell.type in cell.types){
  subset = subset(sc, integrated_annotations == cell.type & condition_collapsed %in% c('mild', 'severe'))
  subset[['ADT_UPDATED']] = NULL
  
  print(DotPlot(subset, features = rev(genes), dot.scale = 8, assay = 'RNA') +
  coord_flip() +
  scale_colour_gradient2(low = 'blue', mid = 'lightgrey', high = 'red', 
                                        midpoint = 0, limits = c(-2,2), oob = scales::squish) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  ggtitle(cell.type))
}
dev.off()



#---- Marker expression

subset = subset(sc, condition %in% c('active_severe', 'active_mild'))
subset$condition = droplevels(subset$condition)

markers = c('CLEC2B', 'KLRF1')
pdf(file = paste0(indir, 'integrated_Tcells_CLEC2B_KLRF1.pdf'))
for (marker in markers){
  print(VlnPlot(subset, feature = marker, split.by = 'condition', 
        pt.size = 0, split = TRUE, cols = viridis(2)) + 
          coord_flip() +
          xlab(''))
}
dev.off()



#---- Heatmap of ADT library expression

DefaultAssay(sc) = 'ADT'
markers = rownames(sc)[24:35]

pdf(file = paste0(indir, 'integrated_Tcells_ADT_selected_expression.pdf'), height = 5.5, width = 8)
DoHeatmap(subset(sc, downsample = 500), 
          features = markers, 
          assay = 'ADT', 
          raster = FALSE, 
          group.colors = cell.type.colors,
          group.bar.height = 0.03, 
          label = FALSE) +
  scale_fill_gradient2(low = viridis(2)[1], 
                       mid = '#F0F0F0', 
                       high = viridis(2)[2], 
                       limits = c(-1,1), 
                       oob = scales::squish)
dev.off()



#---- KLRG1 and IL7R expression 

DefaultAssay(sc) = 'RNA'
genes = c('KLRG1', 'IL7R')
for (gene in genes){
  pdf(file = paste0(indir, 'integrated_Tcells_', gene, '.pdf'), height = 5, width = 5)
  print(VlnPlot(sc, feature = gene, 
          pt.size = 0, 
          group.by = 'integrated_annotations', 
          slot = 'scale.data',
          cols = cell.type.colors,
          sort = TRUE) +
    NoLegend() +
    theme(axis.title.x = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_text(angle = 90)) +
    ylab(paste0(gene, ' expression')) +
    ggtitle(''))
  dev.off()
}



#---- Top5 specificity score (genesorteR) per cluster

sg = run_genesorter(sc, write.file = FALSE)
sg = sg$specScore
genes = list()
n = 5
for (i in 1:ncol(sg)){
  genes[[i]] = rownames(head(sg[order(sg[,i], decreasing = TRUE),], n))
}
top.genes = unique(unlist(genes))


pdf(file = paste0(indir, 'integrated_Tcells_genesorter_top5.pdf'), height = 12, width = 5.5)
DotPlot(sc, 
        features = rev(top.genes), 
        assay = 'RNA') + 
  coord_flip() +
  theme(axis.text.x = element_text(angle = 90, size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank()) +
scale_colour_gradient2(low = viridis(2)[1], mid = '#F5F5F5', high = viridis(2)[2], 
                       midpoint = 0, limits = c(-2,2), oob = scales::squish)
dev.off()
