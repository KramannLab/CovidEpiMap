# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de


#---- Compute p-values for selected PROGENy pathways (active severe vs active mild)

library(dplyr)
library(Seurat)
library(tibble)
library(viridis)
library(tidyverse)
library(tidyr)
source('sc_source/sc_source.R')
indir = '~/sciebo/CovidEpiMap/integrated/'
out.dir =  '~/sciebo/CovidEpiMap/tf_pathway_activity/progeny_sig_test/'
setwd(indir)
sc = readRDS('integrated.RNA.Tcells.annotated.rds')


# Prepare data frame with pathwway scores
progeny_scores_df = as.data.frame(t(GetAssayData(sc, slot = 'scale.data', assay = 'progeny'))) %>%
  rownames_to_column('Cell') %>%
  gather(Pathway, Activity, -Cell)
progeny_scores_df = inner_join(progeny_scores_df, CellsClusters)
progeny_scores_df$condition = sapply(strsplit(progeny_scores_df$CellType,'\\.'), `[`, 1)
progeny_scores_df$cell.type = sapply(strsplit(progeny_scores_df$CellType,'\\.'), `[`, 2)


# Effector memory T cells 1
conditions = c('active_mild', 'active_severe')
celltype = 'CD8+ effector memory T cells 1'
pathways = c('JAK-STAT')
test.stats = compute_stats(df = progeny_scores_df, celltype = celltype, 
                           pathways = pathways, conditions = conditions, plot = TRUE)

file.prefix = paste0(gsub(' ', '_', celltype), '_selected_pvalues_stats')
write.table(test.stats$stats, file = paste0(out.dir, file.prefix, '.txt'),
            sep = '\t',
            row.names = FALSE,
            quote = FALSE)

pdf(file = paste0(out.dir, file.prefix, '.pdf'))
test.stats$p
dev.off()


# Effector memory T cells 2
celltype = 'CD8+ effector memory T cells 2'
pathways = c('JAK-STAT', 'TGFb', 'MAPK')
test.stats = compute_stats(df = progeny_scores_df, celltype = celltype, 
                           pathways = pathways, conditions = conditions, plot = TRUE)

file.prefix = paste0(gsub(' ', '_', celltype), '_selected_pvalues_stats')
write.table(test.stats$stats, file = paste0(out.dir, file.prefix, '.txt'),
            sep = '\t',
            row.names = FALSE,
            quote = FALSE)

pdf(file = paste0(out.dir, file.prefix, '.pdf'))
test.stats$p
dev.off()

# CD8+ TEMRA cells
celltype = 'CD8+ TEMRA cells'
pathways = c('JAK-STAT', 'TNFa', 'NFkB')
test.stats = compute_stats(df = progeny_scores_df, celltype = celltype, 
                           pathways = pathways, conditions = conditions, plot = TRUE)

file.prefix = paste0(gsub(' ', '_', celltype), '_selected_pvalues_stats')
write.table(test.stats$stats, file = paste0(out.dir, file.prefix, '.txt'),
            sep = '\t',
            row.names = FALSE,
            quote = FALSE)

pdf(file = paste0(out.dir, file.prefix, '.pdf'))
test.stats$p
dev.off()


# CD8+ NK−like TEMRA cells
celltype = 'CD8+ NK-like TEMRA cells'
pathways = c('JAK-STAT', 'TNFa', 'NFkB')
test.stats = compute_stats(df = progeny_scores_df, celltype = celltype, 
                           pathways = pathways, conditions = conditions, plot = TRUE)

file.prefix = paste0(gsub(' ', '_', celltype), '_selected_pvalues_stats')
write.table(test.stats$stats, file = paste0(out.dir, file.prefix, '.txt'),
            sep = '\t',
            row.names = FALSE,
            quote = FALSE)

pdf(file = paste0(out.dir, file.prefix, '.pdf'))
test.stats$p
dev.off()


# CD8+ cycling effector T cells
celltype = 'CD8+ cycling effector T cells'
pathways = c('P53', 'JAK-STAT', 'TGFb', 'TNFa', 'NFkB')
test.stats = compute_stats(df = progeny_scores_df, celltype = celltype, 
                           pathways = pathways, conditions = conditions, plot = TRUE)

file.prefix = paste0(gsub(' ', '_', celltype), '_selected_pvalues_stats')
write.table(test.stats$stats, file = paste0(out.dir, file.prefix, '.txt'),
            sep = '\t',
            row.names = FALSE,
            quote = FALSE)

pdf(file = paste0(out.dir, file.prefix, '.pdf'))
test.stats$p
dev.off()


