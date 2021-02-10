# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de


#---- TF activity inference with DoRothEA and pathway activity inference with PROGENy

library(progeny)
library(dorothea)
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
setwd(indir)
sc = readRDS('integrated.RNA.Tcells.annotated.rds')



#---- Run DGEA where we keep all genes and no pAdj cutoff 
sc$celltype.condition = paste(sc$integrated_annotations, sc$condition)
Idents(sc) = 'celltype.condition'


cell.types = names(cell.type.colors)

for (cell.type in cell.types){
	# Healthy vs active mild
	case = 'active_mild'
	control = 'healthy'

	markers = get_markers(sc = sc, condition = case, control = control, 
						cell.type = cell.type, only.sig = FALSE, 
						logfc.threshold = 0, out.dir = 'diff_genes_for_dorothea')


	# Healthy vs active severe
	case = 'active_severe'
	control = 'healthy'

	markers = get_markers(sc = sc, condition = case, control = control, 
						cell.type = cell.type, only.sig = FALSE, 
						logfc.threshold = 0, out.dir = 'diff_genes_for_dorothea')


	# Mild vs severe active
	case = 'active_severe'
	control = 'active_mild'

	markers = get_markers(sc = sc, condition = case, control = control, 
						cell.type = cell.type, only.sig = FALSE, 
						logfc.threshold = 0, out.dir = 'diff_genes_for_dorothea')


	# Mild vs severe recovered
	case = 'recovered_severe'
	control = 'recovered_mild'

	markers = get_markers(sc = sc, condition = case, control = control, 
						cell.type = cell.type, only.sig = FALSE, 
						logfc.threshold = 0, out.dir = 'diff_genes_for_dorothea')
}



#---- Run DoRothEA on contrasts

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



# Healthy vs active mild
diff.indir = '~/sciebo/CovidEpiMap/diff_expression/diff_genes_for_dorothea/'
out.dir =  '~/sciebo/CovidEpiMap/tf_pathway_activity/'
case = 'active_mild'
control = 'healthy'

df = run_dorothea(case = case, control = control, diff.indir = diff.indir, out.dir = out.dir, 
				dorothea_regulon_human = dorothea_regulon_human, regulon = regulon)

pdf(file = paste0(out.dir, case, '_vs_', control, '.pdf'), width = 12, height = 6)
plot_dorothea(df = df, case = case, control = control)
dev.off()


# Healthy vs active severe
case = 'active_severe'
control = 'healthy'

df = run_dorothea(case = case, control = control, diff.indir = diff.indir, out.dir = out.dir, 
                  dorothea_regulon_human = dorothea_regulon_human, regulon = regulon)

pdf(file = paste0(out.dir, case, '_vs_', control, '.pdf'), width = 12, height = 6)
plot_dorothea(df = df, case = case, control = control)
dev.off()


# Mild vs severe active
case = 'active_severe'
control = 'active_mild'

df = run_dorothea(case = case, control = control, diff.indir = diff.indir, out.dir = out.dir, 
                  dorothea_regulon_human = dorothea_regulon_human, regulon = regulon)

pdf(file = paste0(out.dir, case, '_vs_', control, '.pdf'), width = 12, height = 6)
plot_dorothea(df = df, case = case, control = control)
dev.off()


# Mild vs severe recovered
case = 'recovered_severe'
control = 'recovered_mild'

df = run_dorothea(case = case, control = control, diff.indir = diff.indir, out.dir = out.dir, 
                  dorothea_regulon_human = dorothea_regulon_human, regulon = regulon)

pdf(file = paste0(out.dir, case, '_vs_', control, '.pdf'), width = 12, height = 6)
plot_dorothea(df = df, case = case, control = control)
dev.off()



#---- Run PROGENy

sc$integrated_annotations_condition = paste(sc$condition, sc$integrated_annotations, sep = '.')
Idents(sc) = 'integrated_annotations_condition'

sc = progeny(sc, scale = FALSE, organism = 'Human', top = 500, 
			perm = 1, return_assay = TRUE)
sc = ScaleData(sc, assay = 'progeny')


# Create dataframe of clusters
CellsClusters = data.frame(Cell = names(Idents(sc)),
                            CellType = as.character(Idents(sc)),
                            stringsAsFactors = FALSE)

# Transform to data frame
progeny_scores_df = as.data.frame(t(GetAssayData(sc, slot = 'scale.data', assay = 'progeny'))) %>%
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

condition = c('healthy', 'active_mild', 'active_severe',
               'recovered_mild', 'recovered_severe')
annotation$condition = factor(annotation$condition, levels = condition)
condition.colors = viridis(5)
names(condition.colors) = condition

paletteLength = 100
progenyBreaks = c(seq(min(summarized_progeny_scores_df), 0,
                      length.out = ceiling(paletteLength/2) + 1),
                  seq(max(summarized_progeny_scores_df)/paletteLength,
                      max(summarized_progeny_scores_df),
                      length.out = floor(paletteLength/2)))

pdf(file = paste0(out.dir, 'progeny_heatmap.pdf'), width = 10, height = 5)
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

