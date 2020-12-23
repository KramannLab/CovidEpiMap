# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de


library(Seurat)
indir = '~/sciebo/CovidEpiMap/integrated/'
datdir = '~/sciebo/CovidEpiMap/tcr/epitopes/'
source('sc_source/sc_source.R')


sc = readRDS(file = paste0(indir, 'integrated.RNA.Tcells.annotated.rds'))
DefaultAssay(sc) = 'RNA'


# Read epitope-binding/clonotype data
conditions = levels(sc$condition)
tables = lapply(conditions, function(x){
				read.table(file = paste0(datdir, x, '_bycell-selected.csv'), 
					sep = ',', header = TRUE)})


# Subset to only keep cells with enriched dextramer binding compared to negative control
table = do.call(rbind, tables)
table = table[which(table$Log2FC > 0),]