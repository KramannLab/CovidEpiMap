# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de


#---- TF activity inference with DoRothEA and pathway activity inference with PROGENy (A0101-2 binding severe vs mild)

library(progeny)
library(dplyr)
library(Seurat)
library(tibble)
library(viridis)
library(tidyverse)
library(pheatmap)
library(tidyr)
library(viper)
source('sc_source/sc_source.R')
indir = '~/sciebo/CovidEpiMap/integrated/'
outdir = '~/sciebo/CovidEpiMap/epitope_analysis/severe_vs_mild/'


#---- TF activity inference with DoRothEA

# Prepare human DoRothEA regulons
dorothea.path = 'https://raw.githubusercontent.com/saezlab/ConservedFootprints/master/data/dorothea_benchmark/regulons/dorothea_regulon_human_v1.csv'
dorothea_regulon_human = read_csv(dorothea.path)


# Group regulons
regulon = dorothea_regulon_human %>%
  dplyr::filter(confidence %in% c('A','B','C')) %>% 
  split(.$tf) %>%
  map(function(dat) {
    tf = dat %>% distinct(tf) %>% pull()
    targets = setNames(dat$mor, dat$target)
    likelihood = dat$likelihood
    list(tfmode = targets, likelihood = likelihood)
  })

TF_activities_df = data.frame(tf = unique(dorothea_regulon_human$tf), 
                              row.names = unique(dorothea_regulon_human$tf))
cell.types = c('CD8+ TEMRA cells', 'CD8+ effector memory T cells 1')

for (cell.type in cell.types){
  cell.type.name = gsub(' ', '_', cell.type)
  dge = read.table(file = paste0(outdir, cell.type.name, '_A0101-2_binding_severe_vs_mild_dge.txt'),
                   header = TRUE, sep = '\t')
  
  # Estimate z-score values for the gene expression signature (GES)
  myStatistics = matrix(dge$avg_log2FC, dimnames = list(dge$gene, 'avg_log2FC'))
  myPvalue = matrix(dge$p_val, dimnames = list(dge$gene, 'p_val'))
  mySignature = (qnorm(myPvalue/2, lower.tail = FALSE) * sign(myStatistics))[, 1]
  mySignature = mySignature[order(mySignature, decreasing = T)]
  
  # Estimate TF activity
  mrs = msviper(ges = mySignature, regulon = regulon, minsize = 4, 
                ges.filter = FALSE, verbose = FALSE)
  TF_activities = data.frame(Regulon = names(mrs$es$nes),
                             Size = mrs$es$size[ names(mrs$es$nes) ], 
                             NES = mrs$es$nes, 
                             p.value = mrs$es$p.value, 
                             FDR = p.adjust(mrs$es$p.value, method = 'fdr'))
  TF_activities = TF_activities[order(TF_activities$FDR),]
  TF_activities_df[[cell.type]] = TF_activities[rownames(TF_activities_df),'NES']
  
  # Write to file
  write.table(TF_activities, 
              file = paste0(outdir, cell.type.name, '_A0101-2_binding_severe_vs_mild_tf_activity.txt'),
              sep = '\t',
              row.names = FALSE,
              quote = FALSE)
}

pdf(file = paste0(outdir, 'A0101-2_binding_severe_vs_mild_tf_activity.pdf'), height = 5)
plot_dorothea(df = TF_activities_df, case = 'severe', control = 'mild')
dev.off()


#---- Pathway activity inference with PROGENy

# Get unique A0101-2-binding cells
dex.subset = readRDS(file = paste0(indir, 'A0101-2.unique.binders.rds'))

# Subset to relevant cell types and conditions
conditions = c('mild', 'severe')
dex.subset = subset(dex.subset, integrated_annotations %in% cell.types & condition_collapsed %in% conditions)
dex.subset$condition_collapsed = droplevels(dex.subset$condition_collapsed)
dex.subset$integrated_annotations_condition = paste(dex.subset$condition_collapsed, dex.subset$integrated_annotations, sep = '.')
Idents(dex.subset) = 'integrated_annotations_condition'

# Run PROGENy
dex.subset = progeny(dex.subset, scale = FALSE, organism = 'Human', top = 500, 
             perm = 1, return_assay = TRUE)
dex.subset = ScaleData(dex.subset, assay = 'progeny')


# Create dataframe of clusters
CellsClusters = data.frame(Cell = names(Idents(dex.subset)),
                           CellType = as.character(Idents(dex.subset)),
                           stringsAsFactors = FALSE)

# Transform to data frame
progeny_scores_df = as.data.frame(t(GetAssayData(dex.subset, slot = 'scale.data', assay = 'progeny'))) %>%
  rownames_to_column('Cell') %>%
  gather(Pathway, Activity, -Cell)

# Match Progeny scores with the clusters
progeny_scores_df = inner_join(progeny_scores_df, CellsClusters)

# Summarize Progeny scores 
summarized_progeny_scores = progeny_scores_df %>%
  group_by(Pathway, CellType) %>%
  summarise(avg = mean(Activity), std = sd(Activity))

# Create dataframe for plotting
summarized_progeny_scores_df = summarized_progeny_scores %>%
  dplyr::select(-std) %>%
  spread(Pathway, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)


# Plot heatmap
annotation = data.frame(condition = sapply(strsplit(row.names(summarized_progeny_scores_df),'\\.'), `[`, 1),
                        row.names = row.names(summarized_progeny_scores_df))
celltype = sapply(strsplit(row.names(summarized_progeny_scores_df),'\\.'), `[`, 2)

annotation$condition = factor(annotation$condition, levels = conditions)
condition.colors = viridis(length(conditions))
names(condition.colors) = conditions

paletteLength = 100
progenyBreaks = c(seq(min(summarized_progeny_scores_df), 0,
                      length.out = ceiling(paletteLength/2) + 1),
                  seq(max(summarized_progeny_scores_df)/paletteLength,
                      max(summarized_progeny_scores_df),
                      length.out = floor(paletteLength/2)))

pdf(file = paste0(outdir, 'progeny_heatmap.pdf'), width = 4, height = 5)
pheatmap(t(summarized_progeny_scores_df[,-1]),
         fontsize = 9,
         fontsize_row = 9,
         color = colorRampPalette(c('Darkblue', 'white', 'red'))(paletteLength), 
         breaks = progenyBreaks,
         main = '', 
         angle_col = 90,
         treeheight_col = 0,  
         border_color = NA,
         annotation_col = annotation,
         annotation_colors = list(condition = condition.colors),
         labels_col = celltype,
         cluster_cols = FALSE,
         annotation_names_col = FALSE)
dev.off()



#---- Compute p-values for selected PROGENy pathways

# Prepare data frame with pathwway scores
progeny_scores_df = as.data.frame(t(GetAssayData(dex.subset, slot = 'scale.data', assay = 'progeny'))) %>%
  rownames_to_column('Cell') %>%
  gather(Pathway, Activity, -Cell)
progeny_scores_df = inner_join(progeny_scores_df, CellsClusters)
progeny_scores_df$condition = sapply(strsplit(progeny_scores_df$CellType,'\\.'), `[`, 1)
progeny_scores_df$cell.type = sapply(strsplit(progeny_scores_df$CellType,'\\.'), `[`, 2)


# CD8+ effector memory T cells 1
celltype = 'CD8+ effector memory T cells 1'
pathways = c('TNFa', 'JAK-STAT', 'TGFb', 'MAPK')

test.stats = compute_stats(df = progeny_scores_df, celltype = celltype, 
                           pathways = pathways, conditions = conditions, plot = TRUE)

file.prefix = paste0(gsub(' ', '_', celltype), '_selected_pvalues_stats')
write.table(test.stats$stats, file = paste0(outdir, file.prefix, '.txt'),
            sep = '\t',
            row.names = FALSE,
            quote = FALSE)

pdf(file = paste0(outdir, file.prefix, '.pdf'))
test.stats$p
dev.off()


# CD8+ TEMRA cells
celltype = 'CD8+ TEMRA cells'
pathways = c('TNFa', 'JAK-STAT', 'TGFb', 'MAPK')

test.stats = compute_stats(df = progeny_scores_df, celltype = celltype, 
                           pathways = pathways, conditions = conditions, plot = TRUE)

file.prefix = paste0(gsub(' ', '_', celltype), '_selected_pvalues_stats')
write.table(test.stats$stats, file = paste0(outdir, file.prefix, '.txt'),
            sep = '\t',
            row.names = FALSE,
            quote = FALSE)

pdf(file = paste0(outdir, file.prefix, '.pdf'))
test.stats$p
dev.off()




