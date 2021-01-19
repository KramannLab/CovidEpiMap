# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de


#---- Add dextramer-binding data to Seurat object

library(Seurat)
library(dplyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(viridis)
library(tidytext)
indir = '~/sciebo/CovidEpiMap/integrated/'
datdir = '~/sciebo/CovidEpiMap/tcr/'
outdir = '~/sciebo/CovidEpiMap/epitope_analysis/'
'%ni%' = Negate('%in%')
source('sc_source/sc_source.R')


# Format data
sc = readRDS(file = paste0(indir, 'integrated.RNA.Tcells.annotated.rds'))
DefaultAssay(sc) = 'RNA'
Idents(sc) = 'integrated_annotations'


# Add clonotype size
sc$patient_clonotype = str_c(as.character(sc$patient), '_', sc$TCR_clonotype_id)

df = sc@meta.data %>% 
	group_by(patient_clonotype) %>% 
	add_count(name = 'clonotype_size') %>% 
	as.data.frame()

sc$clonotype_size = df$clonotype_size


# Add meta data
conditions = levels(sc$condition)
sc[['COVID']] = ifelse(sc$condition %in% conditions[-1], 'COVID', 'NON-COVID')
sc$condition_collapsed = sub('active_', '', sc$condition)
sc$condition_collapsed = sub('recovered_', '', sc$condition_collapsed)
sc$condition_collapsed = factor(sc$condition_collapsed, levels = c('healthy', 'mild', 'severe'))


# Subset clonotypes to clonotype size > 5 and binding concordance > 30%
clonotype.table = read.table(file = paste0(datdir, 'epitopes_bc0.3.csv'), sep = ',', header = TRUE)
clonotype.table = clonotype.table[clonotype.table$ClonotypeSize > 5,]
clonotype.table = clonotype.table[clonotype.table$BindingConcordance > 0.3,]
dextramers = unique(clonotype.table$Marker)



#---- Add binding information per dextramer

for (dextramer in dextramers){
	# Get clonotypes with enrichment for dextramer
	clonotype = clonotype.table[which(clonotype.table$Marker == dextramer),]$Clonotype

	# Extract cells belonging to clonotypes with enrichment
	df = sc@meta.data
	binding.cells = rownames(df[which(df$patient_clonotype %in% clonotype),])

	# Add meta data
	dextramer = sub('-', '_', dextramer)
	sc[[dextramer]] = ifelse(colnames(sc) %in% binding.cells, 'YES', 'NO')
	df = sc@meta.data


	# Write binding cell overview to table
	# Per cell type
	binding.info.cell.type = df %>% 
				group_by(integrated_annotations) %>% 
				filter(!!sym(dextramer) == 'YES') %>% 
				count(name = dextramer) %>% 
				as.data.frame

	write.table(binding.info.cell.type, 
				file = paste0(outdir, dextramer, '.binding.count.cell.type.txt'),
				sep = '\t',
				row.names = FALSE,
				quote = FALSE)

	# Per cell type, per patient
	binding.info.cell.type.patient = df %>% 
				group_by(patient, integrated_annotations) %>% 
				filter(!!sym(dextramer) == 'YES') %>% 
				count(name = dextramer) %>% 
				as.data.frame

	write.table(binding.info.cell.type.patient, 
				file = paste0(outdir, dextramer, '.binding.count.cell.type.patient.txt'),
				sep = '\t',
				row.names = FALSE,
				quote = FALSE)

	# Per condition
	binding.info.condition = df %>% 
				group_by(condition) %>% 
				filter(!!sym(dextramer) == 'YES') %>% 
				count(name = dextramer) %>% 
				as.data.frame

	write.table(binding.info.condition, 
				file = paste0(outdir, dextramer, '.binding.count.condition.txt'),
				sep = '\t',
				row.names = FALSE,
				quote = FALSE)

	# Per condition, per cell type
	binding.info.condition.cell.type = df %>% 
				group_by(condition, integrated_annotations) %>% 
				filter(!!sym(dextramer) == 'YES') %>% 
				count(name = dextramer) %>% 
				as.data.frame

	write.table(binding.info.condition.cell.type, 
				file = paste0(outdir, dextramer, '.binding.count.condition.cell.type.txt'),
				sep = '\t',
				row.names = FALSE,
				quote = FALSE)

	# Per condition collapsed
	binding.info.condition = df %>% 
				group_by(condition_collapsed) %>% 
				filter(!!sym(dextramer) == 'YES') %>% 
				count(name = dextramer) %>% 
				as.data.frame

	write.table(binding.info.condition, 
				file = paste0(outdir, dextramer, '.binding.count.condition.collapsed.txt'),
				sep = '\t',
				row.names = FALSE,
				quote = FALSE)

	# Per condition collapsed, per cell type
	binding.info.condition.cell.type = df %>% 
				group_by(condition_collapsed, integrated_annotations) %>% 
				filter(!!sym(dextramer) == 'YES') %>% 
				count(name = dextramer) %>% 
				as.data.frame

	write.table(binding.info.condition.cell.type, 
				file = paste0(outdir, dextramer, '.binding.count.condition.collapsed.cell.type.txt'),
				sep = '\t',
				row.names = FALSE,
				quote = FALSE)
}

saveRDS(sc, file = paste0(indir, 'integrated.RNA.Tcells.annotated.rds'))



#---- Count unique binding cells

dextramers = sub('-', '_', dextramers)
dextramers_select = c('A0101_2', 'A0201_4', 'A0201_6')


for (i in 1:length(dextramers_select)){
	dextramer = dextramers_select[i]
	dextramers_remain = dextramers[dextramers %ni% dextramer]
	subset = subset(sc, !!sym(dextramer) == 'YES')

	for (remain_dex in dextramers_remain){
		subset = subset(subset, !!sym(remain_dex) == 'NO' )
	}

	df = subset@meta.data


	# Write binding cell overview to table
	# Per cell type
	binding.info.cell.type = df %>% 
				group_by(integrated_annotations) %>% 
				count(name = dextramer) %>% 
				as.data.frame

	write.table(binding.info.cell.type, 
				file = paste0(outdir, dextramer, '.unique.binding.count.cell.type.txt'),
				sep = '\t',
				row.names = FALSE,
				quote = FALSE)

	# Per cell type, per patient
	binding.info.cell.type.patient = df %>% 
				group_by(patient, integrated_annotations) %>% 
				count(name = dextramer) %>% 
				as.data.frame

	write.table(binding.info.cell.type.patient, 
				file = paste0(outdir, dextramer, '.unique.binding.count.cell.type.patient.txt'),
				sep = '\t',
				row.names = FALSE,
				quote = FALSE)

	# Per condition
	binding.info.condition = df %>% 
				group_by(condition) %>% 
				count(name = dextramer) %>% 
				as.data.frame

	write.table(binding.info.condition, 
				file = paste0(outdir, dextramer, '.unique.binding.count.condition.txt'),
				sep = '\t',
				row.names = FALSE,
				quote = FALSE)

	# Per condition, per cell type
	binding.info.condition.cell.type = df %>% 
				group_by(condition, integrated_annotations) %>% 
				count(name = dextramer) %>% 
				as.data.frame

	write.table(binding.info.condition.cell.type, 
				file = paste0(outdir, dextramer, '.unique.binding.count.condition.cell.type.txt'),
				sep = '\t',
				row.names = FALSE,
				quote = FALSE)

	# Per condition collapsed
	binding.info.condition = df %>% 
				group_by(condition_collapsed) %>% 
				count(name = dextramer) %>% 
				as.data.frame

	write.table(binding.info.condition, 
				file = paste0(outdir, dextramer, '.unique.binding.count.condition.collapsed.txt'),
				sep = '\t',
				row.names = FALSE,
				quote = FALSE)

	# Per condition collapsed, per cell type
	binding.info.condition.cell.type = df %>% 
				group_by(condition_collapsed, integrated_annotations) %>% 
				count(name = dextramer) %>% 
				as.data.frame

	write.table(binding.info.condition.cell.type, 
				file = paste0(outdir, dextramer, '.unique.binding.count.condition.collapsed.cell.type.txt'),
				sep = '\t',
				row.names = FALSE,
				quote = FALSE)
}



#---- Bar chart of unique dextramer binding counts

indir = '~/sciebo/CovidEpiMap/epitope_analysis/binding_counts_unique/'

# Make dextramer unique binding counts data frame
cell.types = names(cell.type.colors)
cell.types = cell.types[cell.types %ni% c('Gamma Delta T cells', 'MAIT cells', 'Atypical NKT cells')]
dextramers = c('A0101_2', 'A0201_4', 'A0201_6')
counts = data.frame(cell.type = rep(cell.types, 3), 
                    condition = rep(c('healthy', 'mild', 'severe'), each = length(cell.types)))


# Add counts from dextramers with unique binding
counts = data.frame()
for (dextramer in dextramers){
  file.prefix = paste0(dextramer, '.unique.binding.count.condition.collapsed.cell.type.txt')
  data = read.table(file = paste0(indir, file.prefix), header = TRUE, sep = '\t')
  data$dextramer = colnames(data)[ncol(data)]
  colnames(data)= c('condition', 'cell.type', 'count', 'dextramer')
  counts = rbind(counts, data)
}
counts = counts[counts$cell.type %in% cell.types,]
counts$condition = factor(counts$condition, levels = c('healthy', 'mild', 'severe'))

# Plot
pdf(file = paste0(indir, 'unique_binding_counts.pdf'), width = 7, height = 4)
ggplot(counts) +
  geom_bar(aes(x = reorder_within(cell.type, -count, dextramer), y = count, fill = condition),
           stat = 'identity') +
  scale_x_reordered() +
  facet_wrap(~ dextramer, scales = 'free') +
  scale_fill_viridis(discrete = TRUE) +
  theme_cowplot() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        axis.text.x = element_text(color = 'black', size = 8, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(color = 'black', size = 8),
        legend.text = element_text(size = 10),
        axis.ticks = element_blank()) +
  ylab('Unique binding count')
dev.off()

