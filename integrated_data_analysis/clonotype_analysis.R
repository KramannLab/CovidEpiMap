# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de


#---- Analysis of TCR clonotype overlap

library(scRepertoire)
library(Seurat)
library(dplyr)
library(stringr)
library(cowplot)
library(vegan)
library(reshape2)
'%ni%' = Negate('%in%')
indir = '~/sciebo/CovidEpiMap/integrated/'
outdir = '~/sciebo/CovidEpiMap/tcr/'
source('sc_source/sc_source.R')

# Get meta data
patients = c('1', '2', '3', 
             '5', '6', '7', 
             '11', '12', '13', 
             '15', '18', '19', 
             '29', '31', '32')
condition = c('active_mild', 'active_mild', 'active_severe', 
              'recovered_mild', 'recovered_mild', 'recovered_mild', 
              'active_severe', 'recovered_severe', 'recovered_severe',
              'active_mild', 'recovered_severe', 'active_severe',
              'healthy', 'healthy', 'healthy')

# Read data
sc = readRDS(file = paste0(indir, 'integrated.RNA.Tcells.annotated.rds'))
DefaultAssay(sc) = 'RNA'
Idents(sc) = 'integrated_annotations'


# Subset to MHC class I/II restricted T cells
sc = subset(sc, integrated_annotations %ni% c('Gamma Delta T cells', 'MAIT cells', 'Atypical NKT cells'))
sc$Clonotype = str_c(sc$patient, '_', sc$TCR_clonotype_id)

sc.2 = RenameCells(sc, 
                   new.names = paste0(sc$orig.ident, '_', 
                                      sc$condition, '_', 
                                      gsub('^.*_GEX_SURF_', '', colnames(sc))))

# Get TCR clonotype information
tcr_list = list()
in.dir = '~/sciebo/CovidEpiMap/data/'
for (patient in patients){
  tcr = read.csv(file = paste0(in.dir, patient, '_TCR/filtered_contig_annotations.csv'))
  tcr_list[[patient]] = tcr
}
combined = combineTCR(tcr_list, samples = paste0(patients, '_GEX_SURF'), ID = condition, cells = 'T-AB')

sce = combineExpression(combined, 
                        sc.2, 
                        cloneCall = 'gene+nt',
                        cloneTypes = c(None = 0, 
                                       Single = 1, 
                                       Small = 5, 
                                       Medium = 20, 
                                       Large = 100, 
                                        Hyperexpanded = 2000))


#  Plot clonal homeostasis
by_cluster = expression2List(sce, group = 'integrated_annotations')
by_condition = expression2List(sce, group = 'condition')
by_cluster_condition = expression2List(sce, group = 'integrated_annotations_condition')

cell.types = levels(droplevels((sc$integrated_annotations)))
levels = paste0(rep(cell.types, each = 5), 
                rep(levels(sc$condition), length(cell.types)))

pdf(file = paste0(outdir, 'clonal_homeostasis_with_patient29.pdf'), height = 5)
clonalHomeostasis(by_condition[levels(sc$condition)], cloneCall = 'gene+nt') +
  scale_fill_manual(name = 'Clonotype Group', values = viridis(5)) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust =1),
        axis.title.x = element_blank(),
        axis.ticks = element_blank())
clonalHomeostasis(by_cluster[cell.types], cloneCall = 'gene+nt') +
  scale_fill_manual(name = 'Clonotype Group', values = viridis(5)) +
  theme_cowplot() +
  theme(axis.text.x = element_text(size = 6, angle = 90, vjust = 1, hjust =1),
        axis.title.x = element_blank(),
        axis.ticks = element_blank())
clonalHomeostasis(by_cluster_condition[levels], cloneCall = 'gene+nt') +
  scale_fill_manual(name = 'Clonotype Group', values = viridis(5)) +
  theme_cowplot() +
  theme(axis.text.x = element_text(size = 6, angle = 90, vjust = 1, hjust =1),
        axis.title.x = element_blank(),
        axis.ticks = element_blank())
dev.off()


# DimPlot with clonal expansion 
slot(sce, 'meta.data')$cloneType = factor(slot(sce, 'meta.data')$cloneType,
                                          levels = c('Hyperexpanded (100 < X <= 2000)', 
                                                     'Large (20 < X <= 100)',
                                                     'Medium (5 < X <= 20)', 
                                                     'Small (1 < X <= 5)',
                                                     'Single (0 < X <= 1)', 
                                                     NA))
pdf(file = paste0(outdir, 'umap_clonotypes_with_patient29.pdf'), width = 10)
DimPlot(sce, group.by = 'cloneType', 
        cols = rev(viridis(5))) +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank())
dev.off()



# Compute morisita overlap index of clonotypes with vegan package
dm = 1 - as.matrix(vegdist(table(sc$integrated_annotations, sc$Clonotype), method = 'horn'))
coef_matrix = melt(get_upper_triangle(dm, diag = FALSE), na.rm = TRUE)

# Plot
pdf(file = paste0(outdir, 'morisita_overlap_index.pdf'), height = 5, width = 6)
upper_trig_tile_plot(coef_matrix) +
  labs(fill = 'Morisita Horn similarity index')
dev.off()

