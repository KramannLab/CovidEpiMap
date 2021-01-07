# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de


#---- Violin plots of canonical cell type markers

library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(viridis)
source('sc_source/sc_source.R')


indir = '~/sciebo/CovidEpiMap/integrated/'
sc = readRDS(file = paste0(indir, 'integrated.RNA.Tcells.annotated.rds'))
DefaultAssay(sc) = 'RNA'
Idents(sc) = 'integrated_annotations'

markers = c(
  'CD8+ naive T cells' = c('CD3', 'CD8A', 'CD8B', 'CD45RA', 'CCR7', 
                          'CD197', 'SELL', 'ACTN1', 'CD27', 'IL7R', 
                          'CD127', 'TCF7', 'CD28'),
  'CD8+ central memory T cells' = c('CD3', 'CD8A', 'CD8B', 'CD45R0', 'CCR7', 
                                    'CD197', 'SELL', 'CD27', 'IL7R', 'CD127', 
                                    'TCF7', 'CD28'),
  'CD8+ CD73+ regulatory T cells' = c('CD3', 'CD8A', 'CD8B', 'CCR7', 'CD27', 
                                      'IL7R', 'CD127', 'TCF7', 'NT5E', 'CCR9'),
  'CD8+ TEMRA cells' = c('CD3', 'CD8A', 'CD8B', 'CD45RA', 'GZMB', 
                        'GZMH', 'GNLY', 'FCGR3A', 'FGFBP2','CX3CR1', 'HLA-DR'),
  'CD8+ NK-like TEMRA cells' = c('CD3', 'CD8A', 'CD8B', 'CD45RA', 'CX3CR1', 
                                'KLRC2', 'KLRF1', 'KIR2DL3', 'KIR3DL2', 'KIR3DL1', 
                                'NCR1', 'CD160'),
  'CD8+ effector memory T cells 1' = c('CD3', 'CD8A', 'CD8B', 'CD45R0', 'CD27', 
                                      'CD28', 'SELL', 'HLA-DR'),
  'CD8+ effector memory T cells 2' = c('CD3', 'CD8A', 'CD8B', 'CD45R0', 'CD27', 
                                      'CD28', 'SELL', 'HLA-DR', 'CD38', 'CX3CR1', 'FGFBP2'),
  'CD8+ cycling effector T cells' = c('CD3', 'CD8A', 'CD8B', 'SELL', 'CD27', 
                                      'CX3CR1', 'HLA-DR', 'MKI67', 'PCNA', 'MCM2'),
  'CD8+ NK-like early effector T cells' = c('CD3', 'CD8A', 'CD8B', 'SELL', 'ACTN1', 
                                            'CD27', 'IL7R', 'CD127', 'TCF7', 'ZNF683', 
                                            'LEF1', 'KLRC2', 'NCR3'),
  'Atypical NKT cells' = c('CD3', 'CD8A', 'CD8B', 'NKG7', 'KLRB1', 
                          'CD161', 'TRAV12-3', 'TRBV7-8', 'TRBV6-5'),
  'CD8+ exhausted T cells' = c('CD3', 'CD8A', 'CD8B', 'CD45RA', 'SELL', 
                              'CD27', 'CD28', 'HLA-DR', 'CD38', 'CTLA4', 
                              'TIGIT', 'CD279', 'HAVCR2', 'PRDM1', 'IL10', 'CXCR6'),
  'MAIT cells' = 'TRAV1-2',
  'Gamma Delta T cells' = 'TRDC')

pdf(file = paste0(indir, 'cell_type_markers.pdf'), width = 14, height = 6)
get_violin(object = sc, features.use = unique(markers))
dev.off()



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
ggplot(df, aes(fill=cluster, y=mean, x=condition)) + 
    geom_bar(stat="identity", position = position_fill(reverse = TRUE)) + 
labs(y = 'Average proportion', x = element_blank(), fill = 'Cell type') +
theme_classic() +
theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
scale_fill_manual(values = cell.type.colors)
dev.off()

pdf(file = paste0(indir, 'integrated_Tcells_barchart_celltype_condition_collapsed_average.pdf'))
ggplot(df, aes(fill=cluster, y=mean, x=condition_collapsed)) + 
    geom_bar(stat="identity", position = position_fill(reverse = TRUE)) + 
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


# Dextramer A0101-2 binding cells

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



#---- Plot selected GO terms

# From DGEA between cell types (active severe vs active mild)
in.dir = '~/sciebo/CovidEpiMap/geneset_enrichment/'
comp = 'active_severe_vs_active_mild'
outdir = paste0(in.dir, comp, '/')

# Define selected GO terms
term.list = list(
  'CD8+ naive T cells' = c('T cell proliferation',
                           'Positive regulation of T cell proliferation',
                           'Positive regulation of viral life cycle',
                           'Positive regulation of immune response'),
  'CD8+ NK-like early effector T cells' = c('T cell receptor signaling pathway',
                                            'Positive regulation of cell activation',
                                            'Regulation of T cell activation',
                                            'Positive regulation of immune response',
                                            'Lymphocyte activation',
                                            'Activation of immune response',
                                            'Innate immune response',
                                            'Adaptive Immune response',
                                            'Cell Killing',
                                            'Antigen processing and presentation of peptide or polysaccharide antigen via MHC class II'),
  'CD8+ effector memory T cells 1' = c('T cell receptor signaling pathway',
                                       'T cell activation',
                                       'Positive regulation of Immune response',
                                       'Lymphocyte activation',
                                       'Antigen receptor mediated signaling pathway',
                                       'Antigen processing and presentation of peptide or polysaccharide antigen via MHC class II',
                                       'Response to type I interferon',
                                       'Negative regulation of viral process',
                                       'Respiratory electron transport chain',
                                       'Oxidative phosphorylation',
                                       'Cellular respiration'),
  'CD8+ effector memory T cells 2' = c('Response to virus',
                                       'Negative regulation of viral process',
                                       'Negative regulation of viral life cycle',
                                       'Response to type I interferon',
                                       'Cellular ketone metabolic process'),
  'CD8+ cycling effector T cells' = c('Viral gene expression',
                                      'Response to type I interferon',
                                      'Positive regulation of canonical Wnt signaling pathway',
                                      'Mitotic cell cycle',
                                      'Cellular ketone metabolic process',
                                      'Cell division',
                                      'Cell cycle process',
                                      'Cell cycle phase transition',
                                      'Cell cycle G2 M Phase transition'),
  'CD8+ NK-like TEMRA cells' = c('Viral genome replication',
                                 'Negative regulation of viral genome replication',
                                 'Negative regulation of viral process',
                                 'Negative regulation of viral life cycle',
                                 'T cell receptor signaling pathway',
                                 'Response to virus',
                                 'Defense response to virus',
                                 'Response to type I interferon',
                                 'Response to interferon beta',
                                 'Antigen processing and presentation of peptide or polysaccharide antigen via MHC class II'), 
  'CD8+ TEMRA cells'= c('Viral genome replication',
                        'Negative regulation of viral genome replication',
                        'Negative regulation of viral process',
                        'Negative regulation of viral life cycle',
                        'T cell activation',
                        'Response to type I interferon',
                        'Regulation of T cell activation',
                        'Positive regulation of immune response',
                        'Antigen processing and presentation of peptide or polysaccharide antigen via MHC class II',
                        'Adaptive immune response'))


width = 7.5
for (cell.type in names(term.list)){
  # Get selected terms
  terms = term.list[[cell.type]]
  
  # Read data
  indir = paste0(in.dir, comp, '/', gsub(' ', '_', cell.type), '/')
  gse = read.table(file = paste0(indir, comp, '_', gsub(' ', '_', cell.type), '_gsea_C5_BP.txt'),
                   header = TRUE, sep = '\t')
  
  # Plot
  if ('Antigen processing and presentation of peptide or polysaccharide antigen via MHC class II' %in% terms){width = 11}
  pdf(file = paste0(outdir, comp, '_', gsub(' ', '_', cell.type), '_selected_GO.pdf'), width = width, height = 3)
  print(plot_go_nice(gse = gse, terms = terms))
  dev.off()
  width = 7.5
}



# From temporal DGEA with TradeSeq
# Start vs end test
in.dir = '~/sciebo/CovidEpiMap/trajectory_analysis/'
comp = 'start_vs_end_test'
outdir = paste0(in.dir, comp, '/')

term.list = list('lineage1' = c('Viral gene expression',
                                'Response to cytokine',
                                'Protein targeting to membrane',
                                'Immune effector process',
                                'Defense response',
                                'Cell activation'),
                 'lineage2' = c('Viral gene expression',
                                'Negative regulation of viral genome replication',
                                'Regulation of Natural Killer cell activation',
                                'Positive regulation of natural killer cell mediated immunity',
                                'Positive regulation of cytokine production involved in immune response',
                                'Positive regulation of innate immune response',
                                'Negative regulation of lymphocyte activation')) 

width = 7.5
for (lineage in names(term.list)){
  # Get selected terms
  terms = term.list[[lineage]]
  
  # Read data
  gse = read.table(file = paste0(in.dir, comp, '/', comp, '_', lineage, '_gsea_C5_BP.txt'),
                   header = TRUE, sep = '\t')
  
  # Plot
  if ('Positive regulation of cytokine production involved in immune response' %in% terms){width = 9}
  pdf(file = paste0(outdir, comp, '_', lineage, '_selected_GO.pdf'), width = width, height = 3)
  print(plot_go_nice(gse = gse, terms = terms))
  dev.off()
  width = 7.5
}


# Lineage split test and diff end test
term.list = list('lineage_split_test' = c('Sensory Perception of Smell',
                                          'Regulation of Natural Killer Cell Activation',
                                          'Positive regulation of natural killer cell mediated immunity',
                                          'Positive regulation of cytokine production involved in immune response',
                                          'Positive regulation of immune effector process',
                                          'Apoptotic cell clearance'),
                 'diff_end_test' = c('Positive regulation of natural killer cell mediated immunity',
                                     'ARP2 3 Complex mediated Actin Nucleation',
                                     'Actin Polymerization or Depolymerization',
                                     'Apoptotic mitochondrial changes',
                                     'Apoptotic process'))

for (comp in names(term.list)){
  outdir = paste0(in.dir, comp, '/')
  
  # Get selected terms
  terms = term.list[[comp]]
  
  # Read data
  gse = read.table(file = paste0(in.dir, comp, '/', comp, '_gsea_C5_BP.txt'),
                   header = TRUE, sep = '\t')
  
  # Plot
  pdf(file = paste0(outdir, comp, '_selected_GO.pdf'), width = 9, height = 3)
  print(plot_go_nice(gse = gse, terms = terms))
  dev.off()
}


# Condition test
comp = 'condition_test'
outdir = paste0(in.dir, comp, '/')

term.list = list('lineage1' = c('Viral gene expression',
                                'Translational initiation',
                                'Cytoplasmatic translation',
                                'Innate Immune response',
                                'Adaptive Immune response',
                                'Immune effector process',
                                'Regulation of immune response'),
                 'lineage2' = c('Viral gene expression',
                                'Regulation of Natural Killer cell activation',
                                'Regulation of cytokine production involved in immune response',
                                'Positive regulation of Natural Killer cell mediated immunity',
                                'Positive regulation of innate immune response',
                                'Positive regulation of cytokine production involved in immune response')) 

width = 7.5
for (lineage in names(term.list)){
  # Get selected terms
  terms = term.list[[lineage]]
  
  # Read data
  gse = read.table(file = paste0(in.dir, comp, '/', comp, '_', lineage, '_gsea_C5_BP.txt'),
                   header = TRUE, sep = '\t')
  
  # Plot
  if ('Positive regulation of cytokine production involved in immune response' %in% terms){width = 9}
  pdf(file = paste0(outdir, comp, '_', lineage, '_selected_GO.pdf'), width = width, height = 3)
  print(plot_go_nice(gse = gse, terms = terms))
  dev.off()
  width = 7.5
}



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

