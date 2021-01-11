# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de


#---- Plot selected GO terms

library(dplyr)
source('sc_source/sc_source.R')


# From DGEA between cell types (active severe vs active mild)
in.dir = '~/sciebo/CovidEpiMap/geneset_enrichment/'
comp = 'active_severe_vs_active_mild'
outdir = paste0(in.dir, comp, '/')

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
  'CD8+ cycling effector T cells' = c('Positive regulation of canonical Wnt signaling pathway',
                                      'Mitotic cell cycle',
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



# From DGEA of epitope A0101-2-binding cells (severe vs mild)
indir = '~/sciebo/CovidEpiMap/epitope_analysis/severe_vs_mild/'
comp = 'severe_vs_mild'

term.list = list('CD8+ effector memory T cells 1' = c('T cell receptor signaling pathway',
                                                    'Cytokine mediated signaling pathway',
                                                    'Activation of Immune response',
                                                    'PID NFAT TFPathway',
                                                    'PID AP1 Pathway'),
                 'CD8+ TEMRA cells' = c('T cell receptor signaling pathway',
                                        'Activation of Immune response',
                                        'PID AP1 Pathway'))

for (cell.type in names(term.list)){
  # Get selected terms
  terms = term.list[[cell.type]]
  
  # Read data
  file.prefix = paste0(gsub(' ', '_', cell.type), '_A0101-2_binding_', comp)
  gse1 = read.table(file = paste0(indir, file.prefix, '_gsea_C5_BP.txt'), header = TRUE, sep = '\t')
  gse2 = read.table(file = paste0(indir, file.prefix, '_gsea_C2_PID.txt'), header = TRUE, sep = '\t')
  gse = rbind(gse1, gse2)
  
  # Plot
  pdf(file = paste0(indir, file.prefix, '_selected_GO.pdf'), width = 6, height = 3)
  print(plot_go_nice(gse = gse, terms = terms))
  dev.off()
}



# From temporal DGEA with TradeSeq
# Start vs end test
in.dir = '~/sciebo/CovidEpiMap/trajectory_analysis/'
comp = 'start_vs_end_test'
outdir = paste0(in.dir, comp, '/')

term.list = list('lineage1' = c('Protein targeting to membrane',
                                'Immune effector process',
                                'Defense response',
                                'Cell activation'),
                 'lineage2' = c('Negative regulation of viral genome replication',
                                'Regulation of Natural Killer cell activation',
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
                                          'Positive regulation of immune effector process',
                                          'Apoptotic cell clearance'),
                 'diff_end_test' = c('ARP2 3 Complex mediated Actin Nucleation',
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

term.list = list('lineage1' = c('Cytoplasmatic translation',
                                'Innate Immune response',
                                'Adaptive Immune response',
                                'Immune effector process',
                                'Regulation of immune response'),
                 'lineage2' = c('Regulation of Natural Killer cell activation',
                                'Regulation of cytokine production involved in immune response',
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

