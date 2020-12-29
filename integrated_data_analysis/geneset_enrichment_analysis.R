# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de


#---- Geneset enrichment analysis of DE genes (active severe vs active mild)

library(Seurat)
indir = '~/sciebo/CovidEpiMap/diff_expression/diff_genes_for_dorothea/'
out.dir = '~/sciebo/CovidEpiMap/geneset_enrichment/'
source('sc_source/sc_source.R')
set.seed(42)

sc = readRDS(file = '~/sciebo/CovidEpiMap/integrated/integrated.RNA.Tcells.annotated.rds')
bg.genes = rownames(sc)

pattern = '.active_severe.vs.active_mild.txt'
de.files = list.files(path = indir, pattern = pattern)

for (de.file in de.files){
  # Read table
  markers = read.table(file = paste0(indir, de.file), 
                     sep = '\t',
                     header = TRUE)
  markers = markers[order(markers$avg_log2FC, decreasing = TRUE),]
  rownames(markers) = markers$gene
  
  # Get cell type
  cell.type = sub(pattern, '', de.file)
  cell.type = sub('integrated.diff.genes.', '', cell.type)
  outdir = paste0(out.dir, 'active_severe_vs_active_mild/', cell.type, '/')
  dir.create(file.path(outdir))
  
  # GSEA
  stats = markers$avg_log2FC
  names(stats) = rownames(markers)
  stats = stats[!is.na(stats)]
  file.prefix = paste0('active_severe_vs_active_mild_', cell.type) 
  
  # GO
  run_gsea(bg.genes = bg.genes, stats = stats, 
           category = 'C5', subcategory = 'BP',
           out.dir = outdir, plot.title = 'GO',
           file.prefix = file.prefix)
  
  # PID
  run_gsea(bg.genes = bg.genes, stats = stats, 
           category = 'C2', subcategory = 'PID',
           out.dir = outdir, plot.title = 'PID',
           file.prefix = file.prefix)
  
  # Immunological signature
  run_gsea(bg.genes = bg.genes, stats = stats, 
           category = 'C7',
           out.dir = outdir, plot.title = 'Immunological Signature',
           file.prefix = file.prefix)
}



#---- Geneset enrichment analysis of DE genes (COVID vs NON-COVID)

# Read data
indir = '~/sciebo/CovidEpiMap/diff_expression/'
markers = read.table(file = paste0(indir, 'integrated.diff.genes.COVID.vs.NON-COVID.txt'), 
                     sep = '\t',
                     header = TRUE)
markers = markers[order(markers$avg_log2FC, decreasing = TRUE),]
rownames(markers) = markers$gene

file.prefix = 'COVID_vs_NON-COVID'
outdir = paste0(out.dir, file.prefix, '/')
dir.create(file.path(outdir))


# GSEA
stats = markers$avg_log2FC
names(stats) = rownames(markers)
stats = stats[!is.na(stats)]

# GO
run_gsea(bg.genes = bg.genes, stats = stats, 
         category = 'C5', subcategory = 'BP',
         out.dir = outdir, plot.title = 'GO',
         file.prefix = file.prefix)

# PID
run_gsea(bg.genes = bg.genes, stats = stats, 
         category = 'C2', subcategory = 'PID',
         out.dir = outdir, plot.title = 'PID',
         file.prefix = file.prefix)

# Immunological signature
run_gsea(bg.genes = bg.genes, stats = stats, 
         category = 'C7',
         out.dir = outdir, plot.title = 'Immunological Signature',
         file.prefix = file.prefix)

