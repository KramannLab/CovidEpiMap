# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de


#---- TF activity inference with DoRothEA and pathway activity inference with PROGENy

library(dorothea)
library(dplyr)
library(Seurat)
library(tibble)
library(pheatmap)
library(tidyverse)
library(tidyr)
library(viper)


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
diff.indir = '~/sciebo/CovidEpiMap/integrated/diff_genes_for_dorothea/'
out.dir =  '~/sciebo/CovidEpiMap/tf_activity_contrast/'
case = 'active_mild'
control = 'healthy'

df = run_dorothea(case = case, control = control, diff.indir = diff.indir, 
				dorothea_regulon_human = dorothea_regulon_human, regulon = regulon)

pdf(file = paste0(out.dir, control, '_vs_',case, '.pdf'), width = 12, height = 6)
plot_dorothea(df = df, case = case, control = control)
dev.off()


# Healthy vs active severe
case = 'active_severe'
control = 'healthy'

df = run_dorothea(case = case, control = control, diff.indir = diff.indir, 
				dorothea_regulon_human = dorothea_regulon_human, regulon = regulon)

pdf(file = paste0(out.dir, control, '_vs_',case, '.pdf'), width = 12, height = 6)
plot_dorothea(df = df, case = case, control = control)
dev.off()


# Mild vs severe active
case = 'active_severe'
control = 'active_mild'

df = run_dorothea(case = case, control = control, diff.indir = diff.indir, 
				dorothea_regulon_human = dorothea_regulon_human, regulon = regulon)

pdf(file = paste0(out.dir, control, '_vs_',case, '.pdf'), width = 12, height = 6)
plot_dorothea(df = df, case = case, control = control)
dev.off()


# Mild vs severe recovered
case = 'recovered_severe'
control = 'recovered_mild'

df = run_dorothea(case = case, control = control, diff.indir = diff.indir, 
				dorothea_regulon_human = dorothea_regulon_human, regulon = regulon)

pdf(file = paste0(out.dir, control, '_vs_',case, '.pdf'), width = 12, height = 6)
plot_dorothea(df = df, case = case, control = control)
dev.off()


     




