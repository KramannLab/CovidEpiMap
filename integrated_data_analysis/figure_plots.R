# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de


#---- Violin plots of canonical cell type markers

library(Seurat)
library(ggplot2)
library(cowplot)
source('sc_source/sc_source.R')


indir = '~/sciebo/CovidEpiMap/integrated/'
sc = readRDS(file = paste0(indir, 'integrated.RNA.Tcells.annotated.rds'))
DefaultAssay(sc) = 'RNA'
Idents(sc) = 'integrated_annotations'



# Combined marker violin plot
get_violin = function(object, features.use){
    p = VlnPlot(object = object, 
                 features = features.use, 
                 pt.size = 0, cols = cell.type.colors,
                 combine = FALSE)
    
    ps = lapply(p, function(x) x + coord_flip() + NoLegend() +
                     theme_bw() +
                     theme(plot.title = element_text(angle = 90), legend.position = 'none',
                           axis.title.x = element_blank(),
                           axis.text.x = element_blank(),
                           axis.ticks.x = element_blank(),
                           axis.text.y = element_blank(),
                           axis.title.y = element_blank(),
                           axis.ticks.y = element_blank(),
                           plot.margin = unit(c(0, 0, 0, 0), 'cm'),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           panel.background = element_rect(colour = 'black', size = 1)))
    
    p = plot_grid(plotlist = ps, nrow = 1, align = 'h')
    
    return(p)
}



# Define markers (mix of ADT and GEX library)
markers = c(
  'CD8+ naive T cells' = c('CD3', 'CD8A', 'CD8B', 
                          'CD45RA', 'CCR7', 'CD197', 
                          'SELL', 'ACTN1', 'CD27', 
                          'IL7R', 'CD127', 'TCF7', 
                          'CD28'),
  'CD8+ central memory T cells' = c('CD3', 'CD8A', 'CD8B',
                                    'CD45R0', 'CCR7', 'CD197', 
                                    'SELL', 'CD27', 'IL7R',
                                    'CD127', 'TCF7', 'CD28'),
  'CD8+ CD73+ regulatory T cells' = c('CD3', 'CD8A', 'CD8B',
                                      'CCR7', 'CD27', 'IL7R',
                                      'CD127', 'TCF7', 'NT5E',
                                      'CCR9'),
  'CD8+ TEMRA cells' = c('CD3', 'CD8A', 'CD8B',
                        'CD45RA', 'GZMB', 'GZMH',
                        'GNLY', 'FCGR3A', 'FGFBP2',
                        'CX3CR1', 'HLA-DR'),
  'CD8+ NK-like TEMRA cells' = c('CD3', 'CD8A', 'CD8B',
                                'CD45RA', 'CX3CR1', 'KLRC2',
                                'KLRF1', 'KIR2DL3', 'KIR3DL2',
                                'KIR3DL1', 'NCR1', 'CD160'),
  'CD8+ effector memory T cells 1' = c('CD3', 'CD8A', 'CD8B',
                                      'CD45R0', 'CD27', 'CD28',
                                      'SELL', 'HLA-DR'),
  'CD8+ effector memory T cells 2' = c('CD3', 'CD8A', 'CD8B',
                                      'CD45R0', 'CD27', 'CD28',
                                      'SELL', 'HLA-DR', 'CD38',
                                      'CX3CR1', 'FGFBP2'),
  'CD8+ cycling effector T cells' = c('CD3', 'CD8A', 'CD8B',
                                      'SELL', 'CD27', 'CX3CR1',
                                      'HLA-DR', 'MKI67', 'PCNA', 
                                      'MCM2'),
  'CD8+ NK-like early effector T cells' = c('CD3', 'CD8A', 'CD8B',
                                            'SELL', 'ACTN1', 'CD27',
                                            'IL7R', 'CD127', 'TCF7', 
                                            'ZNF683', 'LEF1', 'KLRC2', 
                                            'NCR3'),
  'Atypical NKT cells' = c('CD3', 'CD8A', 'CD8B',
                          'NKG7', 'KLRB1', 'CD161',
                          'TRAV12-3', 'TRBV7-8', 'TRBV6-5'),
  'CD8+ exhausted T cells' = c('CD3', 'CD8A', 'CD8B',
                              'CD45RA', 'SELL', 'CD27',
                              'CD28', 'HLA-DR', 'CD38',
                              'CTLA4', 'TIGIT', 'CD279',
                              'HAVCR2', 'PRDM1', 'IL10',
                              'CXCR6'),
  'MAIT cells' = 'TRAV1-2',
  'Gamma Delta T cells' = 'TRDC')


p = get_violin(object = sc, features.use = unique(markers))

pdf(file = paste0(indir, 'cell_type_markers.pdf'), width = 14, height = 6)
p
dev.off()



#---- DotPlot of top10 genes based on specificity score per cluster

# Genesorter
sg = run_genesorter(sc, assay = 'RNA', slot = 'data', write.file = TRUE,
  out.dir = indir, file.name.prefix = 'integrated.tcells.')
specScore = sg$specScore


# Dotplot of top 10 genes per cluster
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



