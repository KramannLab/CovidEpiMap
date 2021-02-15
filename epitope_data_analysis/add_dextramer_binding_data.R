# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de


#---- Add dextramer-binding data to Seurat object

# Source for binding concordance (10x Genomics)
# https://www.10xgenomics.com/resources/application-notes/a-new-way-of-exploring-immunity-linking-highly-multiplexed-antigen-recognition-to-immune-repertoire-and-phenotype/

library(Seurat)
library(dplyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(viridis)
library(tidytext)
library(reshape2)
library(matrixStats)
indir = '~/sciebo/CovidEpiMap/integrated/'
outdir = '~/sciebo/CovidEpiMap/epitope_analysis/'
'%ni%' = Negate('%in%')
source('sc_source/sc_source.R')


# Format data
sc.all = readRDS(file = paste0(indir, 'integrated.RNA.Tcells.annotated.rds'))
DefaultAssay(sc.all) = 'RNA'
Idents(sc.all) = 'integrated_annotations'

# Subset to relevant cell types
excluded.cell.types = c('Gamma Delta T cells', 'MAIT cells', 'Atypical NKT cells')
sc = subset(sc.all, integrated_annotations %ni% excluded.cell.types)
sc$integrated_annotations = droplevels(sc$integrated_annotations)
sc$Clonotype = str_c(sc$patient, '_', sc$TCR_clonotype_id)

# Add meta data
conditions = levels(sc$condition)
sc$condition_collapsed = sub('active_', '', sc$condition)
sc$condition_collapsed = sub('recovered_', '', sc$condition_collapsed)
sc$condition_collapsed = factor(sc$condition_collapsed, levels = c('healthy', 'mild', 'severe'))


# Get dextramer counts and convert to log2 scale (add pseudocount)
adt = GetAssayData(sc, assay = 'ADT_UPDATED', slot = 'counts')
adt_log2 = as.matrix(log2(adt + 1))

# Median expression of negative control
negative_controls = c('A0101-33', 'A0201-34', 'A0301-35', 'general-36')
negctrl = as.array(colMedians(adt_log2[negative_controls,]))

# L2FC between the dextramers (n = 15) and the median of the negative controls
df_scores = sweep(adt_log2[1:15,], 2, negctrl, '-')

# Get clonotype size
clonotype_size = table(sc$Clonotype)

# Add clonotype and size to L2FCs
df_scores = melt(df_scores, na.rm = T)
colnames(df_scores) = c('dextramer', 'barcode', 'L2FC')
df_scores$clonotype = sc$Clonotype[as.character(df_scores$barcode)]
df_scores$clonotype_size = as.numeric(clonotype_size[df_scores$clonotype])

# Define binding concordance as percentage of cells within a clonotype
# with a L2FC > 2 between the dextramer and the median of negative control
# Keep clonotypes with binding concordance > 0.3 and size > 5
binding_concordance = df_scores %>%
  group_by(clonotype, dextramer) %>%
  mutate(binding_concordance = sum(L2FC > 2)/clonotype_size) %>%
  filter(binding_concordance > 0.3,
         clonotype_size > 5)

# Get binding clonotypes
binding_concordance = binding_concordance %>% 
  select(-barcode, -L2FC) %>% 
  as.data.frame %>% 
  unique



#---- Add binding information per dextramer

dextramers = as.vector(unique(binding_concordance$dextramer))
outdir = '~/sciebo/CovidEpiMap/epitope_analysis/binding_counts/'

for (dex in dextramers){
  # Get clonotypes with enrichment for dextramer
  clonotype = binding_concordance %>% 
    filter(dextramer == dex) %>% 
    select(clonotype) %>% 
    unique
  clonotype = clonotype$clonotype
  
  # Extract cells belonging to clonotypes with enrichment
  df = sc@meta.data
  binding.cells = rownames(df[which(df$Clonotype %in% clonotype),])
  
  # Add meta data
  dex = sub('-', '_', dex)
  sc[[dex]] = ifelse(colnames(sc) %in% binding.cells, 'YES', 'NO')
  df = sc@meta.data
  
  # Write binding cell overview to table
  # Per cell type, per patient
  binding.info.cell.type.patient = df %>% 
    group_by(patient, integrated_annotations) %>% 
    filter(!!sym(dex) == 'YES') %>% 
    dplyr::count(name = dex) %>% 
    as.data.frame
  
  write.table(binding.info.cell.type.patient, 
              file = paste0(outdir, dex, '.binding.count.cell.type.patient.txt'),
              sep = '\t',
              row.names = FALSE,
              quote = FALSE)
  
  # Per condition collapsed
  binding.info.condition = df %>% 
    group_by(condition_collapsed) %>% 
    filter(!!sym(dex) == 'YES') %>% 
    dplyr::count(name = dex) %>% 
    as.data.frame
  
  write.table(binding.info.condition, 
              file = paste0(outdir, dex, '.binding.count.condition.collapsed.txt'),
              sep = '\t',
              row.names = FALSE,
              quote = FALSE)
  
  # Per condition collapsed, per cell type
  binding.info.condition.cell.type = df %>% 
    group_by(condition_collapsed, integrated_annotations) %>% 
    filter(!!sym(dex) == 'YES') %>% 
    dplyr::count(name = dex) %>% 
    as.data.frame
  
  write.table(binding.info.condition.cell.type, 
              file = paste0(outdir, dex, '.binding.count.condition.collapsed.cell.type.txt'),
              sep = '\t',
              row.names = FALSE,
              quote = FALSE)
}


# Save data
dextramers = sub('-', '_', dextramers)
df = sc@meta.data
df.all = sc.all@meta.data

for (dex in dextramers){
  df.sub = df[,dex]
  names(df.sub) = rownames(df)
  sc.all[[dex]] = df.sub[rownames(df.all)]
}

saveRDS(sc.all, file = paste0(indir, 'integrated.RNA.Tcells.annotated.rds'))



#---- Count unique binding cells

dextramers_select = c('A0101_2', 'A0201_4', 'A0201_6', 'A1101_29')
outdir = '~/sciebo/CovidEpiMap/epitope_analysis/binding_counts_unique/'

for (i in 1:length(dextramers_select)){
  dextramer = dextramers_select[i]
  dextramers_remain = dextramers[dextramers %ni% dextramer]
  subset = subset(sc, !!sym(dextramer) == 'YES')
  
  for (remain_dex in dextramers_remain){
    subset = subset(subset, !!sym(remain_dex) == 'NO' )
  }
  
  df = subset@meta.data
  
  # Write binding cell overview to table
  # Per cell type, per patient
  binding.info.cell.type.patient = df %>% 
    group_by(patient, integrated_annotations) %>% 
    dplyr::count(name = dextramer) %>% 
    as.data.frame
  
  write.table(binding.info.cell.type.patient, 
              file = paste0(outdir, dextramer, '.unique.binding.count.cell.type.patient.txt'),
              sep = '\t',
              row.names = FALSE,
              quote = FALSE)
  
  # Per condition collapsed
  binding.info.condition = df %>% 
    group_by(condition_collapsed) %>% 
    dplyr::count(name = dextramer) %>% 
    as.data.frame
  
  write.table(binding.info.condition, 
              file = paste0(outdir, dextramer, '.unique.binding.count.condition.collapsed.txt'),
              sep = '\t',
              row.names = FALSE,
              quote = FALSE)
  
  # Per condition collapsed, per cell type
  binding.info.condition.cell.type = df %>% 
    group_by(condition_collapsed, integrated_annotations) %>% 
    dplyr::count(name = dextramer) %>% 
    as.data.frame
  
  write.table(binding.info.condition.cell.type, 
              file = paste0(outdir, dextramer, '.unique.binding.count.condition.collapsed.cell.type.txt'),
              sep = '\t',
              row.names = FALSE,
              quote = FALSE)
}
