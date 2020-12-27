# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de


#---- TF activity inference with DoRothEA and pathway activity inference with PROGENy

#devtools::install_github('saezlab/progeny')
library(progeny)
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
out.dir =  '~/sciebo/CovidEpiMap/tf_pathway_activity/'
case = 'active_mild'
control = 'healthy'

df = run_dorothea(case = case, control = control, diff.indir = diff.indir, 
				dorothea_regulon_human = dorothea_regulon_human, regulon = regulon)

pdf(file = paste0(out.dir, control, '_vs_',case, '.pdf'), width = 12, height = 6)
plot_dorothea(df = df, case = case, control = control)
dev.off()

# Write to file
df = df[complete.cases(df),]
write.table(df, file = paste0(out.dir, control, '_vs_',case, '.txt'),
			sep = '\t',
			row.names = FALSE,
			quote = FALSE)


# Healthy vs active severe
case = 'active_severe'
control = 'healthy'

df = run_dorothea(case = case, control = control, diff.indir = diff.indir, 
				dorothea_regulon_human = dorothea_regulon_human, regulon = regulon)

pdf(file = paste0(out.dir, control, '_vs_',case, '.pdf'), width = 12, height = 6)
plot_dorothea(df = df, case = case, control = control)
dev.off()

# Write to file
df = df[complete.cases(df),]
write.table(df, file = paste0(out.dir, control, '_vs_',case, '.txt'),
			sep = '\t',
			row.names = FALSE,
			quote = FALSE)


# Mild vs severe active
case = 'active_severe'
control = 'active_mild'

df = run_dorothea(case = case, control = control, diff.indir = diff.indir, 
				dorothea_regulon_human = dorothea_regulon_human, regulon = regulon)

pdf(file = paste0(out.dir, control, '_vs_',case, '.pdf'), width = 12, height = 6)
plot_dorothea(df = df, case = case, control = control)
dev.off()

# Write to file
df = df[complete.cases(df),]
write.table(df, file = paste0(out.dir, control, '_vs_',case, '.txt'),
			sep = '\t',
			row.names = FALSE,
			quote = FALSE)


# Mild vs severe recovered
case = 'recovered_severe'
control = 'recovered_mild'

df = run_dorothea(case = case, control = control, diff.indir = diff.indir, 
				dorothea_regulon_human = dorothea_regulon_human, regulon = regulon)

pdf(file = paste0(out.dir, control, '_vs_',case, '.pdf'), width = 12, height = 6)
plot_dorothea(df = df, case = case, control = control)
dev.off()

# Write to file
df = df[complete.cases(df),]
write.table(df, file = paste0(out.dir, control, '_vs_',case, '.txt'),
			sep = '\t',
			row.names = FALSE,
			quote = FALSE)
     


#---- Run PROGENy  cell clusters

cell.clusters = data.frame(cell = names(Idents(sc)), 
							cell.type = as.character(Idents(sc)),
							stringsAsFactors = FALSE)

sc = progeny(sc, scale = FALSE, organism = 'Human', top = 500, 
			perm = 1, return_assay = TRUE)

sc = ScaleData(sc, assay = 'progeny')

progeny.scores = as.data.frame(t(GetAssayData(sc, slot = 'scale.data',
					assay = 'progeny'))) %>%
				rownames_to_column('cell') %>%
				gather(Pathway, Activity, -cell)

progeny.scores = inner_join(progeny.scores, cell.clusters)

summarized.progeny.scores = progeny.scores %>%
							group_by(Pathway, cell.type) %>%
							summarise(avg = mean(Activity), std = sd(Activity))

summarized.progeny.df = summarized.progeny.scores %>%
						dplyr::select(-std) %>%
						spread(Pathway, avg) %>%
						data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)


pdf(file = paste0(outdir,'progeny_violin.pdf'), width = 17)
for (gene in progeny.genes){
	print(VlnPlot(sc, features = gene, 
		group.by = 'integrated_annotations', 
		pt.size = 0, split.by = 'condition',
		assay = 'progeny', slot = 'scale.data'))
}
dev.off()



