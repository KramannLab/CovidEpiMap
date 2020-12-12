# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de


#---- Differential gene expression analysis across conditions of integrated data

library(Seurat)
library(dplyr)
library(ggplot2)
indir = '~/Dropbox/CovidEpiMap/integrated/'
setwd(indir)
source('../sc_source/sc_source.R')


sc = readRDS(file = 'integrated.RNA.Tcells.annotated.rds')
DefaultAssay(sc) = 'RNA'

cell.types = levels(Idents(sc))
sc$celltype.condition = paste(Idents(sc), sc$condition)
Idents(sc) = 'celltype.condition'



for (cell.type in cell.types){
	# Healthy vs active mild
	condition = 'active_mild'
	control = 'healthy'

	markers = get_markers(sc, condition = condition, 
		control = control, cell.type = cell.type)

	top_heatmap(sc, markers = markers, cell.type = cell.type, 
		disease = condition, control = control, n = 15)


	# Healthy vs active severe
	condition = 'active_severe'
	control = 'healthy'

	markers = get_markers(sc, condition = condition, 
		control = control, cell.type = cell.type)
	
	top_heatmap(sc, markers = markers, cell.type = cell.type, 
		disease = condition, control = control, n = 15)


	# Mild vs severe active
	condition = 'active_severe'
	control = 'active_mild'

	markers = get_markers(sc, condition = condition, 
		control = control, cell.type = cell.type)
	
	top_heatmap(sc, markers = markers, cell.type = cell.type, 
		disease = condition, control = control, n = 15)


	# Mild vs severe recovered
	condition = 'recovered_severe'
	control = 'recovered_mild'

	markers = get_markers(sc, condition = condition, 
		control = control, cell.type = cell.type)
	
	top_heatmap(sc, markers = markers, cell.type = cell.type, 
		disease = condition, control = control, n = 15)
}


