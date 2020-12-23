# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de


library(Seurat)
indir = '~/sciebo/CovidEpiMap/integrated/'
datdir = '~/sciebo/CovidEpiMap/tcr/epitopes/'
source('sc_source/sc_source.R')


sc = readRDS(file = paste0(indir, 'integrated.RNA.Tcells.annotated.rds'))
DefaultAssay(sc) = 'RNA'
Idents(sc) = 'integrated_annotations'
conditions = levels(sc$condition)


# Read epitope-binding data (single-cell-level)
cell.tables = lapply(conditions, function(x){
				read.table(file = paste0(datdir, x, '_bycell-selected.csv'), 
					sep = ',', header = TRUE)})
cell.table = do.call(rbind, cell.tables)


# Read epitope-binding data (clonotype-level)
# Only clonotypes with size > 5 are included in the tables
clonotype.tables = lapply(conditions, function(x){
				read.table(file = paste0(datdir, x, '_byclone-selected.csv'), 
					sep = ',', header = TRUE)})
clonotype.table = do.call(rbind, clonotype.tables)



# Subset to only keep cell with enriched dextramer-binding compared to negative control (L2FC > 2)

table = table[which(table$Log2FC > 2),]


