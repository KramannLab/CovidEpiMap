# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de


#---- Differential gene expression analysis across conditions of integrated data

library(Seurat)
library(dplyr)
library(ggplot2)
indir = '~/sciebo/CovidEpiMap/integrated/'
outdir = '~/sciebo/CovidEpiMap/diff_expression/diff_genes/'
source('sc_source/sc_source.R')


sc = readRDS(file = paste0(indir, 'integrated.RNA.Tcells.annotated.rds'))
DefaultAssay(sc) = 'RNA'


cell.types = levels(Idents(sc))
sc$celltype.condition = paste(Idents(sc), sc$condition)
Idents(sc) = 'celltype.condition'


# Cell type vs cell type across conditions
for (cell.type in cell.types){
	# Healthy vs active mild
	condition = 'active_mild'
	control = 'healthy'

	markers = get_markers(sc, condition = condition, 
		control = control, 
		cell.type = cell.type,
		out.dir = outdir)

	# Healthy vs active severe
	condition = 'active_severe'
	control = 'healthy'

	markers = get_markers(sc, condition = condition, 
	                      control = control, 
	                      cell.type = cell.type,
	                      out.dir = outdir)

	# Mild vs severe active
	condition = 'active_severe'
	control = 'active_mild'

	markers = get_markers(sc, condition = condition, 
	                      control = control, 
	                      cell.type = cell.type,
	                      out.dir = outdir)

	# Mild vs severe recovered
	condition = 'recovered_severe'
	control = 'recovered_mild'
	
	markers = get_markers(sc, condition = condition, 
	                      control = control, 
	                      cell.type = cell.type,
	                      out.dir = outdir)
}

