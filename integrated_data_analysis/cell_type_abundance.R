# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de


#---- Cell-type abundance testing

# Source
# https://www.nxn.se/valent/2020/11/28/s9jjv32ogiplagwx8xrkjk532p7k28
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggrepel)
library(emmeans)
indir = '~/sciebo/CovidEpiMap/integrated/'
outdir = '~/sciebo/CovidEpiMap/cell_type_abundance/'
source('sc_source/sc_source.R')


sc = readRDS(file = paste0(indir, 'integrated.RNA.Tcells.annotated.rds'))
DefaultAssay(sc) = 'RNA'


# Count total number of cells per patient and number of cells per cell type
df = sc@meta.data
data = df %>%
select(condition, patient, integrated_annotations) %>%
group_by(patient) %>%
count(integrated_annotations, name = 'count', .drop = FALSE) %>%
add_count(wt = count, name = 'total') %>%
mutate(other = total - count) %>%
as.data.frame()


# Add condition to table again
condition = rep(c('healthy', 'active_mild',
				'active_severe', 'recovered_mild',
				'recovered_severe'), 
				each = 3*13)
data = cbind(data, condition)


# Save counts
write.table(data, file = paste0(outdir, 'integrated.cell.type.counts.txt'),
			sep = '\t', row.names = FALSE, quote = FALSE)



# Healthy vs active mild
conditions = c('healthy', 'active_mild')
plots = get_abundance(df = data, condition.vector = conditions)

pdf(file = paste0(outdir, conditions[1], '_vs_', conditions[2], '.pdf'))
plots$p1
plots$p2
dev.off()


# Healthy vs active severe
conditions = c('healthy', 'active_severe')
plots = get_abundance(df = data, condition.vector = conditions)

pdf(file = paste0(outdir, conditions[1], '_vs_', conditions[2], '.pdf'))
plots$p1
plots$p2
dev.off()


# Healthy vs recovered mild
conditions = c('healthy', 'recovered_mild')
plots = get_abundance(df = data, condition.vector = conditions)

pdf(file = paste0(outdir, conditions[1], '_vs_', conditions[2], '.pdf'))
plots$p1
plots$p2
dev.off()


# Healthy vs recovered severe
conditions = c('healthy', 'recovered_severe')
plots = get_abundance(df = data, condition.vector = conditions)

pdf(file = paste0(outdir, conditions[1], '_vs_', conditions[2], '.pdf'))
plots$p1
plots$p2
dev.off()


# Active mild vs active severe
conditions = c('active_mild', 'active_severe')
plots = get_abundance(df = data, condition.vector = conditions)

pdf(file = paste0(outdir, conditions[1], '_vs_', conditions[2], '.pdf'))
plots$p1
plots$p2
dev.off()


# Healthy vs mild
data.coll = data 
data.coll$condition = sub('active_', '', data.coll$condition)
data.coll$condition = sub('recovered_', '', data.coll$condition)

conditions = c('healthy', 'mild')
plots = get_abundance(df = data.coll, condition.vector = conditions)

pdf(file = paste0(outdir, conditions[1], '_vs_', conditions[2], '.pdf'))
plots$p1
plots$p2
dev.off()


# Healthy vs severe
conditions = c('healthy', 'severe')
plots = get_abundance(df = data.coll, condition.vector = conditions)

pdf(file = paste0(outdir, conditions[1], '_vs_', conditions[2], '.pdf'))
plots$p1
plots$p2
dev.off()


# Mild vs severe
conditions = c('mild', 'severe')
plots = get_abundance(df = data.coll, condition.vector = conditions)

pdf(file = paste0(outdir, conditions[1], '_vs_', conditions[2], '.pdf'))
plots$p1
plots$p2
dev.off()

