# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de


#---- Plot selected GO terms

library(dplyr)
library(cowplot)
library(ggplot2)
source('sc_source/sc_source.R')


# From DGEA between cell types (active severe vs active mild)
in.dir = '~/sciebo/CovidEpiMap/geneset_enrichment/'
comp = 'active_severe_vs_active_mild'

term.list = list(
  'CD8+ NK-like early effector T cells' = c('GO T cell receptor signaling pathway',
                                            'GO Positive regulation of cell activation',
                                            'GO Regulation of T cell activation',
                                            'GO Positive regulation of immune response',
                                            'GO Lymphocyte activation',
                                            'GO Activation of immune response',
                                            'GO Innate immune response',
                                            'GO Adaptive Immune response',
                                            'GO Cell Killing',
                                            'GO Antigen processing and presentation of peptide or polysaccharide antigen via MHC class II'),
  'CD8+ effector memory T cells 1' = c('GO T cell receptor signaling pathway',
                                       'GO T cell activation',
                                       'GO Positive regulation of Immune response',
                                       'GO Lymphocyte activation',
                                       'GO Antigen receptor mediated signaling pathway',
                                       'GO Antigen processing and presentation of peptide or polysaccharide antigen via MHC class II',
                                       'GO Response to type I interferon',
                                       'GO Negative regulation of viral process',
                                       'GO Respiratory electron transport chain',
                                       'GO Oxidative phosphorylation',
                                       'GO Cellular respiration',
                                       'GO Response to interferon beta'),
  'CD8+ effector memory T cells 2' = c('GO Response to virus',
                                       'GO Negative regulation of viral process',
                                       'GO Negative regulation of viral life cycle',
                                       'GO Response to type I interferon',
                                       'GO Cellular ketone metabolic process',
                                       'GO Response to type I interferon'),
  'CD8+ cycling effector T cells' = c('GO Positive regulation of canonical Wnt signaling pathway',
                                      'GO Mitotic cell cycle',
                                      'GO Cell division',
                                      'GO Cell cycle process',
                                      'GO Cell cycle phase transition',
                                      'GO Cell cycle G2 M Phase transition',
                                      'GO Response to type I interferon'),
  'CD8+ NK-like TEMRA cells' = c('GO Response to virus',
                                 'GO Defense response to virus',
                                 'GO Response to type I interferon',
                                 'GO Response to interferon beta',
                                 'GO Antigen processing and presentation of peptide or polysaccharide antigen via MHC class II'), 
  'CD8+ TEMRA cells'= c('GO Response to type I interferon',
                        'GO Antigen processing and presentation of peptide or polysaccharide antigen via MHC class II',
                        'GO Adaptive immune response'))

# List for selected GOs
plot.list = list()
# List for all sig GOs (for supplements)
all.plot.list = list()

for (cell.type in names(term.list)){
  # Get selected terms
  terms = term.list[[cell.type]]
  
  # Read data
  indir = paste0(in.dir, comp, '/', gsub(' ', '_', cell.type), '/')
  gse = read.table(file = paste0(indir, comp, '_', gsub(' ', '_', cell.type), '_gsea_C5_BP.txt'),
                   header = TRUE, sep = '\t')
  gse = gse[complete.cases(gse),]
  
  title = paste0(comp, '_', gsub(' ', '_', cell.type))
  plot.list[[title]] = plot_go_nice(gse = gse, terms = terms, title = title)
  
  if (nrow(gse[gse$padj < 0.05,]) > 0){
    gse = gse[gse$padj < 0.05,]
    gse = gse[order(gse$NES, decreasing= TRUE),]
    terms = gse$pathway
    all.plot.list[[title]] = plot_go_nice(gse = gse, terms = head(terms,30), title = title)
  }
}


# From DGEA of epitope A0101-2-binding cells (severe vs mild)
indir = '~/sciebo/CovidEpiMap/epitope_analysis/severe_vs_mild/GSEA/'
comp = 'severe_vs_mild'

term.list = list('CD8+ effector memory T cells 1' = c('GO T cell receptor signaling pathway',
                                                      'GO Cytokine mediated signaling pathway',
                                                      'GO Activation of Immune response',
                                                      'GO Response to Interferon Gamma',
                                                      'GO Interferon Gamma mediated signaling pathway',
                                                      'GO Immune response regulating signaling pathway',
                                                      'GO Antigen receptor mediated signaling pathway',
                                                      'PID NFAT TFPathway',
                                                      'PID AP1 Pathway',
                                                      'PID IL6 7 Pathway'),
                 'CD8+ TEMRA cells' = c('GO T cell receptor signaling pathway',
                                        'GO Activation of Immune response',
                                        'GO Response to Interleukin 1',
                                        'GO Positive regulation of Immune response',
                                        'GO Activation of Innate Immune response',
                                        'GO Innate Immune response activating signal transduction',
                                        'GO Antigen receptor mediated signaling pathway',
                                        'PID AP1 Pathway'))

for (cell.type in names(term.list)){
  # Get selected terms
  terms = term.list[[cell.type]]
  
  # Read data
  title = paste0(gsub(' ', '_', cell.type), '_A0101-2_binding_', comp)
  gse1 = read.table(file = paste0(indir, title, '_gsea_C5_BP.txt'), header = TRUE, sep = '\t')
  gse2 = read.table(file = paste0(indir, title, '_gsea_C2_PID.txt'), header = TRUE, sep = '\t')
  gse = rbind(gse1, gse2)
  gse = gse[complete.cases(gse),]
  
  plot.list[[title]] = plot_go_nice(gse = gse, terms = terms, title = title)
  
  if (nrow(gse[gse$padj < 0.05,]) > 0){
    gse = gse[gse$padj < 0.05,]
    gse = gse[order(gse$NES, decreasing= TRUE),]
    terms = gse$pathway
    all.plot.list[[title]] = plot_go_nice(gse = gse, terms = head(terms,30), title = title)
  }
}


# From temporal DGEA with TradeSeq
# Start vs end test
in.dir = '~/sciebo/CovidEpiMap/trajectory_analysis/'
comp = 'start_vs_end_test'

term.list = list('lineage1' = c('GO Viral gene expression',
                                'GO Response to cytokine',
                                'GO Protein targeting to membrane',
                                'GO Immune effector process',
                                'GO Defense response',
                                'GO Cell activation',
                                'GO Apoptotic process'),
                 'lineage2' = c('GO Regulation of Natural Killer cell activation',
                                'GO Positive regulation of natural killer cell mediated immunity',
                                'GO Positive regulation of cytokine production involved in immune response')) 

for (lineage in names(term.list)){
  # Get selected terms
  terms = term.list[[lineage]]
  
  # Read data
  title = paste0(comp, '_', lineage)
  gse = read.table(file = paste0(in.dir, comp, '/', title, '_gsea_C5_BP.txt'),
                   header = TRUE, sep = '\t')
  gse = gse[complete.cases(gse),]
  
  plot.list[[title]] = plot_go_nice(gse = gse, terms = terms, title = title)
  
  if (nrow(gse[gse$padj < 0.05,]) > 0){
    gse = gse[gse$padj < 0.05,]
    gse = gse[order(gse$NES, decreasing= TRUE),]
    terms = gse$pathway
    all.plot.list[[title]] = plot_go_nice(gse = gse, terms = head(terms,30), title = title)
  }
}


# Lineage split test
term.list = list('lineage_split_test' = c('GO Regulation of Natural Killer Cell Activation',
                                          'GO Positive regulation of natural killer cell mediated immunity',
                                          'GO Positive regulation of cytokine production involved in immune response'))

for (comp in names(term.list)){
  # Get selected terms
  terms = term.list[[comp]]
  
  # Read data
  gse = read.table(file = paste0(in.dir, comp, '/', comp, '_gsea_C5_BP.txt'),
                   header = TRUE, sep = '\t')
  gse = gse[complete.cases(gse),]
  
  plot.list[[comp]] = plot_go_nice(gse = gse, terms = terms, title = comp)
  
  if (nrow(gse[gse$padj < 0.05,]) > 0){
    gse = gse[gse$padj < 0.05,]
    gse = gse[order(gse$NES, decreasing= TRUE),]
    terms = gse$pathway
    all.plot.list[[comp]] = plot_go_nice(gse = gse, terms = head(terms,30), title = comp)
  }
}


# Condition test
comp = 'condition_test'

term.list = list('lineage1' = c('GO Viral gene expression',
                                'GO Translational initiation',
                                'GO Cytoplasmic translation',
                                'GO Innate Immune response',
                                'GO Defense response'),
                 'lineage2' = c('GO Regulation of Natural Killer cell activation',
                                'GO Positive regulation of Natural Killer cell mediated immunity',
                                'GO Positive regulation of cytokine production involved in immune response',
                                'GO Regulation of Natural Killer cell mediated immunity'))


for (lineage in names(term.list)){
  # Get selected terms
  terms = term.list[[lineage]]
  
  # Read data
  title = paste0(comp, '_', lineage)
  gse = read.table(file = paste0(in.dir, comp, '/', title, '_gsea_C5_BP.txt'),
                   header = TRUE, sep = '\t')
  gse = gse[complete.cases(gse),]
  
  plot.list[[title]] = plot_go_nice(gse = gse, terms = terms, title = title)
  
  if (nrow(gse[gse$padj < 0.05,]) > 0){
    gse = gse[gse$padj < 0.05,]
    gse = gse[order(gse$NES, decreasing= TRUE),]
    terms = gse$pathway
    all.plot.list[[title]] = plot_go_nice(gse = gse, terms = head(terms,30), title = title)
  }
}


# Plot selected GOs
outdir = '~/sciebo/CovidEpiMap/geneset_enrichment/'
pdf(file = paste0(outdir, 'selected_GOs_all.pdf'), height = 40, width = 8)
plot_grid(plotlist = plot.list, align = 'hv', ncol = 1)
dev.off()

# Plot all significant GOs
pdf(file = paste0(outdir, 'significant_GOs_all.pdf'), height = 100, width = 25)
plot_grid(plotlist = all.plot.list, align = 'hv', ncol = 1)
dev.off()

