# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de


#--- Add cell type annotations to subset for relevant T cells for per patient clonotyping

library(Seurat)
indir = '~/sciebo/CovidEpiMap/per_patient/'
setwd(indir)
'%ni%' = Negate('%in%')


# Patient 1
patient = '1'
sample = paste0(patient, '_GEX_SURF')
sc = readRDS(file = paste0(sample, '/', sample, '.rds'))

bad.clusters = c('8')
sc = subset(sc, subset = seurat_clusters %ni% bad.clusters)

sc = RenameIdents(sc, `0` = 'CD8+ T cells 1', `1` = 'CD8+ T EMRA cells 1',
		`2` = 'MAIT cells', `3` = 'CD8+ T cells 2',
		`4` = 'CD8+ T EMRA cells 2', `5` = 'CD8+ T cells 3',
		`6` = 'Monocytes/Macrophages', `7` = 'CD8+ central memory T cells',
		`9` = 'CD8+ effector memory T cells', `10` = 'CD8+ T EMRA cells 3',
		`11` = 'CD8+ naive T cells', `12` = 'B cells',
		`13` = 'Gamma/Delta T cells')
sc$orig_annotations = Idents(sc)
saveRDS(sc, file = paste0(sample, '/', sample, 'annotated.rds'))



# Patient 2
patient = '2'
sample = paste0(patient, '_GEX_SURF')
sc = readRDS(file = paste0(sample, '/', sample, '.rds'))

bad.clusters = c('7')
sc = subset(sc, subset = seurat_clusters %ni% bad.clusters)

sc = RenameIdents(sc, `0` = 'CD8+ central memory T cells 1', `1` = 'CD8+ central memory T cells 2',
		`2` = 'CD8+ naive T cells', `3` = 'CD8+ central memory T cells 3',
		`4` = 'CD4+/CD8+ T cells', `5` = 'CD8+ central memory T cells 4',
		`6` = 'CD8+ T cells 1', `8` = 'CD8+ T cells 2',
		`9` = 'B cells', `10` = 'CD8+ T cells 3',
		`11` = 'CD8+ T cells 4', `12` = 'CD8+ T cells 5',
		`13` = 'Gamma/Delta T cells', `14` = 'Monocytes/Macrophages')
sc$orig_annotations = Idents(sc)
saveRDS(sc, file = paste0(sample, '/', sample, 'annotated.rds'))



# Patient 3 
patient = '3'
sample = paste0(patient, '_GEX_SURF')
sc = readRDS(file = paste0(sample, '/', sample, '.rds'))

sc = RenameIdents(sc, `0` = 'CD8+ naive T cells', `1` = 'Gamma/Delta T cells',
		`2` = 'CD8+ T EMRA cells', `3` = 'MAIT cells', 
		`4` = 'CD8+ T cells 1', `5`= 'CD8+ T cells 2',
		`6` = 'CD8+ T cells 3', `7` = 'CD8+ T cells 4',
		`8` = 'CD8+ central memory T cells 1', `9` = 'B cells',
		`10` = 'CD8+ central memory T cells 2', `11` = 'Monocytes/Macrophages',
		`12` = 'Platelets', `13` = 'CD4+/CD8+ T cells')
sc$orig_annotations = Idents(sc)
saveRDS(sc, file = paste0(sample, '/', sample, 'annotated.rds'))



# Patient 5 
patient = '5'
sample = paste0(patient, '_GEX_SURF')
sc = readRDS(file = paste0(sample, '/', sample, '.rds'))

sc = RenameIdents(sc, `0` = 'CD8+ naive T cells 1', `1` = 'MAIT cells',
		`2` = 'CD8+ central memory T cells', `3` = 'CD8+ T cells 1',
		`4` = 'B cells 1', `5` = 'CD8+ T cells 2',
		`6` = 'CD8+ T cells 3', `7` = 'B cells 2',
		`8` = 'Monocytes/Macrophages', `9` = 'CD8+ T cells 4',
		`10` = 'CD8+ naive T cells 2')
sc$orig_annotations = Idents(sc)
saveRDS(sc, file = paste0(sample, '/', sample, 'annotated.rds'))



# Patient 6
patient = '6'
sample = paste0(patient, '_GEX_SURF')
sc = readRDS(file = paste0(sample, '/', sample, '.rds'))

sc = RenameIdents(sc, `0` = 'CD8+ T cells 1', `1` = 'CD8+ central memory T cells 1',
		`2` = 'CD8+ naive T cells 1', `3` = 'CD8+ T cells 2',
		`4` = 'CD8+ T cells 3', `5` = 'CD8+ T cells 4',
		`6` = 'CD8+ effector memory T cells 1', `7` = 'MAIT cells',
		`8` = 'B cells 1', `9` = 'CD8+ naive T cells 2',
		`10` = 'B cells 2', `11` = 'CD8+ T cells 5',
		`12` = 'CD8+ effector memory T cells 2', `13` = 'CD8+ central memory T cells 2',
		`14` = 'B cells 3', `15` = 'B cells 4',
		`16` = 'Monocytes/Macrophages', `17` = 'CD8+ T cells 6',
		`18` = 'Gamma/Delta T cells')
sc$orig_annotations = Idents(sc)
saveRDS(sc, file = paste0(sample, '/', sample, 'annotated.rds'))



# Patient 7
patient = '7'
sample = paste0(patient, '_GEX_SURF')
sc = readRDS(file = paste0(sample, '/', sample, '.rds'))

sc = RenameIdents(sc, `0` = 'MAIT cells 1', `1` = 'CD8+ central memory T cells 1',
		`2` = 'CD8+ naive T cells', `3` = 'CD8+ T cells 1',
		`4` = 'B cells 1', `5` = 'CD8+ T cells 2',
		`6` = 'CD8+ T cells 3', `7` = 'CD8+ T cells 4',
		`8` = 'CD8+ T cells 5', `9` = 'CD8+ T cells 6',
		`10` = 'B cells 2', `11` = 'CD8+ T cells 7',
		`12` = 'CD8+ T EMRA cells', `13` = 'CD8+ T cells 8',
		`14` = 'MAIT cells 2', `15` = 'CD8+ T cells 9',
		`16` = 'CD8+ T cells 10', `17` = 'CD8+ T cells 11',
		`18` = 'CD4+ T cells', `19` = 'Monocytes/Macrophages',
		`20` = 'CD8+ central memory T cells 2')
sc$orig_annotations = Idents(sc)
saveRDS(sc, file = paste0(sample, '/', sample, 'annotated.rds'))



# Patient 11
patient = '11'
sample = paste0(patient, '_GEX_SURF')
sc = readRDS(file = paste0(sample, '/', sample, '.rds'))

sc = RenameIdents(sc, `0` = 'CD8+ T cells 1', `1` = 'CD8+ T cells 2',
		`2` = 'CD8+ T cells 3', `3` = 'CD8+ T cells 4',
		`4` = 'CD8+ T cells 5', `5` = 'CD8+ T cells 6',
		`6` = 'CD8+ T cells 7', `7` = 'B cells',
		`8` = 'CD8+ T cells 8', `9` = 'CD8+ naive T cells',
		`10` = 'CD8+ T cells 9', `11` = 'CD8+ T cells 10',
		`12` = 'CD4+/CD8+ T cells', `13` = 'CD8+ T cells 11',
		`14` = 'CD8+ T cells 12', `15` = 'Monocytes/Macrophages',
		`16` = 'CD8+ T cells 13')
sc$orig_annotations = Idents(sc)
saveRDS(sc, file = paste0(sample, '/', sample, 'annotated.rds'))



# Patient 12
patient = '12'
sample = paste0(patient, '_GEX_SURF')
sc = readRDS(file = paste0(sample, '/', sample, '.rds'))

bad.clusters = c('14')
sc = subset(sc, subset = seurat_clusters %ni% bad.clusters)

sc = RenameIdents(sc, `0` = 'Monocytes/Macrophages 1', `1` = 'Monocytes/Macrophages 2',
		`2` = 'Monocytes/Macrophages 3', `3` = 'Monocytes/Macrophages 4',
		`4` = 'Monocytes/Macrophages 5', `5` = 'Monocytes/Macrophages 6',
		`6` = 'CD8+ T EMRA cells 1', `7` = 'CD8+ T EMRA cells 2',
		`8` = 'Monocytes/Macrophages 7', `9` = 'CD8+ T cells 1',
		`10` = 'B cells 1', `11` = 'Monocytes/Macrophages 8',
		`12` = 'Monocytes/Macrophages 9', `13` = 'CD8+ T EMRA cells 3',
		`15` = 'CD8+ effector memory T cells', `16` = 'B cells 2',
		`17` = 'CD8+ T EMRA cells 4', `18` = 'MAIT cells',
		`19` = 'CD8+ T cells 2')
sc$orig_annotations = Idents(sc)
saveRDS(sc, file = paste0(sample, '/', sample, 'annotated.rds'))



# Patient 13
patient = '13'
sample = paste0(patient, '_GEX_SURF')
sc = readRDS(file = paste0(sample, '/', sample, '.rds'))

bad.clusters = c('12', '15')
sc = subset(sc, subset = seurat_clusters %ni% bad.clusters)

sc = RenameIdents(sc, `0` = 'CD8+ T cells 1', `1` = 'CD8+ T cells 2',
		`2` = 'CD8+ T cells 3', `3` = 'CD8+ T cells 4',
		`4` = 'CD8+ T cells 5', `5` = 'CD8+ naive T cells',
		`6` = 'CD8+ T cells 6', `7` = 'CD8+ T cells 7', 
		`8` = 'CD8+ central memory T cells 1', `9` = 'CD8+ T cells 8',
		`10` = 'CD8+ central memory T cells 2', `11` = 'CD8+ T cells 9',
		`13` = 'Monocytes/Macrophages', `14` = 'CD8+ T cells 10')
sc$orig_annotations = Idents(sc)
saveRDS(sc, file = paste0(sample, '/', sample, 'annotated.rds'))



# Patient 15
patient = '15'
sample = paste0(patient, '_GEX_SURF')
sc = readRDS(file = paste0(sample, '/', sample, '.rds'))

bad.clusters = c('15', '17')
sc = subset(sc, subset = seurat_clusters %ni% bad.clusters)

sc = RenameIdents(sc, `0` = 'CD8+ T EMRA cells 1', `1` = 'CD8+ T EMRA cells 2',
		`2` = 'CD8+ T EMRA cells 3', `3` = 'CD8+ T EMRA cells 4',
		`4` = 'CD8+ effector memory T cells 1', `5` = 'CD8+ T EMRA cells 5',
		`6` = 'CD8+ T EMRA cells 6', `7` = 'B cells',
		`8` = 'CD8+ effector memory T cells 2', `9` = 'CD8+ effector memory T cells 3',
		`10` = 'CD8+ effector memory T cells 4', `11` = 'CD8+ T cells 1',
		`12` = 'CD8+ naive T cells', `13` = 'Monocytes/Macrophages',
		`14` = 'CD8+ T cells 2', `16` = 'Gamma/Delta T cells')
sc$orig_annotations = Idents(sc)
saveRDS(sc, file = paste0(sample, '/', sample, 'annotated.rds'))



# Patient 18
patient = '18'
sample = paste0(patient, '_GEX_SURF')
sc = readRDS(file = paste0(sample, '/', sample, '.rds'))

sc = RenameIdents(sc, `0` = 'B cells 1', `1` = 'CD8+ naive T cells',
		`2` = 'CD8+ T cells 1', `3` = 'B cells 2',
		`4` = 'CD8+ T cells 2', `5` = 'CD8+ T cells 3',
		`6` = 'B cells 3', `7` = 'CD8+ T cells 4',
		`8` = 'CD8+ central memory T cells 1', `9` = 'CD8+ central memory T cells 2',
		`10` = 'B cells 4', `11` = 'MAIT cells', 
		`12` = 'CD8+ T cells 5', `13` = 'CD8+ T cells 6',
		`14` = 'CD8+ effector memory T cells 1', `15` = 'Monocytes/Macrophages',
		`16` = 'Gamma/Delta T cells', `17` = 'CD8+ T cells 7',
		`18` = 'CD8+ effector memory T cells 2')
sc$orig_annotations = Idents(sc)
saveRDS(sc, file = paste0(sample, '/', sample, 'annotated.rds'))



# Patient 19
patient = '19'
sample = paste0(patient, '_GEX_SURF')
sc = readRDS(file = paste0(sample, '/', sample, '.rds'))

sc = RenameIdents(sc, `0` = 'CD8+ T cells 1', `1` = 'CD8+ T cells 2',
		`2` = 'CD8+ T cells 3', `3` = 'CD8+ T EMRA cells 1',
		`4` = 'CD8+ T EMRA cells 2', `5` = 'B cells 1',
		`6` = 'CD8+ effector memory T cells 1', `7` = 'CD8+ effector memory T cells 2',
		`8` = 'CD8+ T cells 4', `9` = 'CD8+ T cells 5',
		`10` = 'CD8+ T EMRA cells 3', `11` = 'CD8+ T EMRA cells 4',
		`12` = 'CD8+ T EMRA cells 5', `13` = 'CD8+ T cells 6',
		`14` = 'CD8+ T cells 7', `15` = 'CD8+ T EMRA cells 6',
		`16` = 'CD8+ T EMRA cells 7', `17` = 'CD8+ effector memory T cells 3',
		`18` = 'CD8+ effector memory T cells 4', `19` = 'CD8+ T EMRA cells 8',
		`20` = 'Monocytes/Macrophages', `21` = 'CD8+ effector memory T cells 5',
		`22` = 'CD8+ effector memory T cells 6', `23` = 'B cells 2',
		`24` = 'CD8+ T cells 8', `25` = 'CD8+ T cells 9',
		`26` = 'CD8+ T EMRA cells 9')
sc$orig_annotations = Idents(sc)
saveRDS(sc, file = paste0(sample, '/', sample, 'annotated.rds'))



# Patient 31
patient = '31'
sample = paste0(patient, '_GEX_SURF')
sc = readRDS(file = paste0(sample, '/', sample, '.rds'))

sc = RenameIdents(sc, `0` = 'B cells 1', `1` = 'CD8+ T cells 1',
		`2` = 'B cells 2', `3` = 'B cells 3',
		`4` = 'CD8+ T EMRA cells 1', `5` = 'CD8+ effector memory T cells 1',
		`6` = 'CD8+ T EMRA cells 2', `7` = 'CD8+ T cells 2',
		`8` = 'B cells 4', `9` = 'CD8+ T cells 3',
		`10` = 'CD8+ T cells 4', `11` = 'B cells 5',
		`12` = 'CD8+ effector memory T cells 2', `13` = 'CD8+ T cells 5',
		`14` = 'B cells 6', `15` = 'CD8+ naive T cells',
		`16` = 'B cells 7', `17` = 'CD8+ effector memory T cells 3',
		`18` = 'CD8+ effector memory T cells 4', `19` = 'CD4+ T cells',
		`20` = 'MAIT cells', `21` = 'B cells 8',
		`22` = 'Monocytes/Macrophages', `23` = 'B cells 9',
		`24` = 'B cells 10', `25` = 'Gamma/Delta T cells',
		`26` = 'B cells 11', `27` = 'B cells 12',
		`28` = 'B cells 13', `29` =  'CD141+ myeloid dendritic cells',
		`30` = 'Non-classical Monocytes/Macrophages')
sc$orig_annotations = Idents(sc)
saveRDS(sc, file = paste0(sample, '/', sample, 'annotated.rds'))



# Patient 32
patient = '32'
sample = paste0(patient, '_GEX_SURF')
sc = readRDS(file = paste0(sample, '/', sample, '.rds'))

sc = RenameIdents(sc, `0` = 'MAIT cells', `1` = 'CD8+ central memory T cells',
		`2` = 'CD8+ T EMRA cells 1', `3` = 'CD8+ T cells 1',
		`4` = 'B cells', `5` = 'CD8+ T EMRA cells 2',
		`6` = 'CD8+ T cells 2', `7` = 'Gamma/Delta T cells')
sc$orig_annotations = Idents(sc)
saveRDS(sc, file = paste0(sample, '/', sample, 'annotated.rds'))


