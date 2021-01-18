# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de


#---- Plot TF activity (active severe vs active mild)

library(progeny)
library(dorothea)
library(dplyr)
library(tibble)
library(viridis)
library(tidyverse)
library(pheatmap)
library(tidyr)
library(viper)
'%ni%' = Negate('%in%')
indir = '~/sciebo/CovidEpiMap/tf_pathway_activity/active_severe_vs_active_mild/'


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

case = 'active_severe'
control = 'active_mild'
diff.indir = '~/sciebo/CovidEpiMap/diff_expression/diff_genes_for_dorothea/'
out.dir =  '~/sciebo/CovidEpiMap/tf_pathway_activity/'

tf.df = run_dorothea(case = case, control = control, diff.indir = diff.indir, out.dir = out.dir, 
                  dorothea_regulon_human = dorothea_regulon_human, regulon = regulon)


# Get significant TFs
pattern = 'active_severe_vs_active_mild_tf_activity_'
case.files = list.files(indir, pattern = pattern)

tf.tables = list()
for (case in case.files){
  cell.type = sub(pattern, '', case)
  cell.type = sub('.txt', '', cell.type)

  tf.table = read.table(file = paste0(indir, case), header = TRUE, sep = '\t')
  tf.table$cell.type = cell.type
  tf.tables[[cell.type]] = tf.table
}
tf.tables = do.call(rbind, tf.tables)
rownames(tf.tables) = NULL

# Highlight the significant TFs
sig.regulons = tf.tables %>% 
  filter(FDR < 0.05) %>% 
  select(Regulon) %>% 
  unique 
sig.regulons = sig.regulons$Regulon

# Get cell types with more than one sig TF
sig.cell.types = tf.tables %>% 
  filter(FDR < 0.05) %>% 
  count(cell.type) %>% 
  filter(n > 1) %>% 
  select(cell.type)
sig.cell.types = sig.cell.types$cell.type


# Format data for plotting
df = tf.df
df$tf = NULL
df = apply(df, 2, as.numeric)
rownames(df) = rownames(tf.df)
df.plot = df[,colnames(df) %in% sig.cell.types]
# Keep only complete cases except if TF is significant in at least one case
df.plot = df.plot[rowSums(!is.na(df.plot)) == ncol(df.plot) | rownames(df.plot) %in% sig.regulons, ]

# Define colour palette
paletteLength = 100
myColor = colorRampPalette(c(viridis(paletteLength)[1], 
                             'white', 
                             viridis(paletteLength)[paletteLength]))(paletteLength)
viperBreaks = c(seq(min(df.plot, na.rm = TRUE), 0, 
                    length.out = ceiling(paletteLength/2) + 1),
                seq(max(df.plot, na.rm = TRUE)/paletteLength, 
                    max(df.plot, na.rm = TRUE), 
                    length.out = floor(paletteLength/2)))

# Cluster complete cases of TFs (otherwise hclust will return error)
df.plot.sub = df.plot[complete.cases(df.plot),]
celltype_cluster = hclust(dist(t(df.plot.sub)))
celltype_ordered = celltype_cluster$labels[celltype_cluster$order]
tf_cluster = hclust(dist(df.plot.sub))
tf_ordered = tf_cluster$labels[tf_cluster$order]

# Manually edit order of top TFs to made output look nicer
top = c('RFX5','RFXANK','RFXAP','STAT1','STAT2','IRF2','JUNB','STAT3','MYC','ZEB2')
tf_ordered = tf_ordered[tf_ordered %ni% top]
# Append TFs with NAs to order of the complete TF
tf_ordered_all = c(top, tf_ordered, rownames(df.plot)[!rownames(df.plot) %in% c(top, tf_ordered)])

# Plot
pdf(file = paste0(indir,'active_severe_active_mild_tf_activity.pdf'), width = 2.3, height = 6)
pheatmap(df.plot[tf_ordered_all,celltype_ordered],
         fontsize_col = 7, 
         fontsize_row = 2.3,
         color = myColor, 
         breaks = viperBreaks, 
         main = '',
         angle_col = 90,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         treeheight_row = 0, 
         treeheight_col = 0,
         labels_row = ifelse(tf_ordered_all %in% sig.regulons, tf_ordered_all, ''))
dev.off()


