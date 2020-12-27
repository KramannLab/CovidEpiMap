# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de


#---- Analysis of dextramer-epitope binding cells 

library(Seurat)
library(dplyr)
library(msigdbr)
library(fgsea)
indir = '~/sciebo/CovidEpiMap/integrated/'
datdir = '~/sciebo/CovidEpiMap/tcr/'
outdir = '~/sciebo/CovidEpiMap/epitope_analysis/'
source('sc_source/sc_source.R')



# GSEA function
run_gsea = function(object, markers, category, plot.title = NULL, subcategory = NULL, out.dir = '.', file.prefix){
	# Fetch geneset
	geneSets = msigdbr(species = 'Homo sapiens', category = category, subcategory = subcategory)
	geneSets = geneSets[geneSets$human_gene_symbol %in% rownames(object),]
	m_list = geneSets %>% split(x = .$human_gene_symbol, f = .$gs_name)

	# Run GSEA
	stats = markers$avg_log2FC
	names(stats) = rownames(markers)
	gsea = fgsea(pathways = m_list, stats = stats, minSize = 10, eps = 0.0)
	order = order(gsea$padj, decreasing = FALSE)

	# Plot
	file.name = paste0(out.dir, '/', file.prefix, '_gsea_', paste0(c(category, subcategory), collapse = '_'))

	pdf(file = paste0(file.name, '.pdf'), width = 10, height = 9)
	print(plot_go(gsea.res = gsea,
			gsea.res.order = order, 
			n = 30, 
			plot.title = plot.title))
	dev.off()

	# Write to file
	write.table(gsea[order, -8], 
			file = paste0(file.name, '.txt'), 
			sep = '\t', 
			row.names = FALSE, 
			quote = FALSE)
}



# Format data
sc = readRDS(file = paste0(indir, 'integrated.RNA.Tcells.annotated.rds'))
DefaultAssay(sc) = 'RNA'
Idents(sc) = 'integrated_annotations'
conditions = levels(sc$condition)

sc$patient_clonotype = paste0(sc$patient, '_', sc$TCR_clonotype_id)
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



#---- Plot dextramer binding on UMAP

dextramers = sub('-', '_', dextramers)

pdf(file = paste0(outdir, 'integrated_Tcells_dextramer_binding.pdf'))
for (dextramer in dextramers){
	print(DimPlot(sc, group.by = dextramer, cols = rev(viridis(2)), order = TRUE) +
		ggtitle(dextramer))
}
dev.off()


# Plot unique A0101-2 binding
dex.subset = subset(sc, A0201_5 == 'NO' & A1101_30 == 'NO' & A0201_6 == 'NO' & 
					A0101_2 == 'YES' & A1101_29 == 'NO' & A1101_23 == 'NO' & 
					A0201_4 == 'NO' & A0201_12 == 'NO')

pdf(file = paste0(outdir, 'integrated_Tcells_dextramer_A0101-2_unique_binding.pdf'))
DimPlot(sc, cells.highlight = colnames(dex.subset),
		cols.highlight = viridis(1), sizes.highlight = 0.1)
dev.off()



#---- Plot clonal expansion of non-binding cells

df = sc@meta.data %>% 
	group_by(patient_clonotype) %>% 
	add_count(name = 'clonotype_size') %>% 
	as.data.frame()

sc$clonotype_size = df$clonotype_size
subset = subset(sc, A0201_5 == 'NO' & A1101_30 == 'NO' & A0201_6 == 'NO' & 
					A0101_2 == 'NO' & A1101_29 == 'NO' & A1101_23 == 'NO' & 
					A0201_4 == 'NO' & A0201_12 == 'NO')

pdf(file = paste0(outdir, 'clonal_expansion_non_binding_cells.pdf'))
FeaturePlot(subset, feature = 'clonotype_size', cols = c('yellow', 'red'), order = TRUE)
dev.off()



#---- Plot clonal expansion in binding cells


pdf(file = paste0(outdir, 'clonal_expansion_binding_cells_covid_split.pdf'), width = 10)
FeaturePlot(subset, feature = 'clonotype_size', cols = c('yellow', 'red'), 
			split.by = 'COVID', max.cutoff = 500, order = TRUE)
dev.off()



#---- DGEA/GSEA of binding COVID vs NON-COVID cells

Idents(subset) = 'COVID'

markers = FindMarkers(subset, ident.1 = 'COVID',
      				ident.2 = 'NON-COVID', 
      				min.pct = 0.25,
      				logfc.threshold = 0)

markers = markers[order(markers$avg_log2FC, decreasing = TRUE),]

# Write to file
markers.out = markers[markers$p_val_adj < 0.05,]
markers.out$gene = rownames(markers.out)

write.table(markers.out[,c(6,1:5)], 
        file = paste0(outdir, 'COVID_vs_NON-COVID_binding_cells_dge.txt'),
        quote = FALSE, 
        sep = '\t', 
        row.names = FALSE)


# GSEA
# GO
run_gsea(object = subset, markers = markers, category = 'C5', subcategory = 'BP',
		out.dir = outdir, plot.title = 'GO',
		file.prefix = paste0('COVID_vs_NON-COVID_binding_cells'))

# PID
run_gsea(object = subset, markers = markers, category = 'C2', subcategory = 'PID',
		out.dir = outdir, plot.title = 'PID',
		file.prefix = paste0('COVID_vs_NON-COVID_binding_cells'))

# Immunological signature
run_gsea(object = subset, markers = markers, category = 'C7',
		out.dir = outdir, plot.title = 'Immunological Signature',
		file.prefix = paste0('COVID_vs_NON-COVID_binding_cells'))




#---- A0101_2-binding vs non-binding cells (active severe)

dextramer = 'A0101_2'
df = sc@meta.data

# Only evaluate cell types with >200 cells in each group
cell.types = c('CD8+ TEMRA cells', 'CD8+ effector memory T cells 1',
				'CD8+ effector memory T cells 2')

# Active severe binding CD8+ TEMRA cells: 2216
# Active severe non-binding CD8+ TEMRA cells: 1280

# Active severe binding CD8+ effector memory T cells 1: 601
# Active severe non-binding CD8+ effector memory T cells 1: 1008

# Active severe binding CD8+ effector memory T cells 2: 417
# Active severe non-binding CD8+ effector memory T cells 2: 247


for (cell.type in cell.types){
	# DGEA
	subset = subset(sc, condition == 'active_severe' & integrated_annotations == cell.type)
	Idents(subset) = dextramer

	markers = FindMarkers(subset, ident.1 = 'YES', ident.2 = 'NO', logfc.threshold = 0)
	markers = markers[order(markers$avg_log2FC, decreasing = TRUE),]


	# Write to file
	markers.out = markers[markers$p_val_adj < 0.05,]
	markers.out$gene = rownames(markers.out)
	cell.type.name = gsub(' ', '_', cell.type)

	write.table(markers.out[,c(6, 1:5)], 
			file = paste0(outdir, dextramer, '/active_severe_', cell.type.name, '_binding_vs_nonbinding_dge.txt'), 
			sep = '\t',
			row.names = FALSE, 
			quote = FALSE)


	# GSEA
	# GO
	run_gsea(object = sc, markers = markers, category = 'C5', subcategory = 'BP',
			out.dir = paste0(outdir, dextramer), cell.type = cell.type, 
			file.prefix = paste0('active_severe_', cell.type.name, '_binding_vs_nonbinding'))

	# PID
	run_gsea(object = sc, markers = markers, category = 'C2', subcategory = 'PID',
			out.dir = paste0(outdir, dextramer), cell.type = cell.type, 
			file.prefix = paste0('active_severe_', cell.type.name, '_binding_vs_nonbinding'))

	# Immunological signature
	run_gsea(object = sc, markers = markers, category = 'C7',
			out.dir = paste0(outdir, dextramer), cell.type = cell.type,
			file.prefix = paste0('active_severe_', cell.type.name, '_binding_vs_nonbinding'))	
}



#---- A0101_2-binding recovered severe vs recovered mild CD8+ TEMRA cells 

# Recovered mild binding CD8+ TEMRA cells: 599
# Recovered severe binding CD8+ TEMRA cells: 250


# DGEA
cell.type = 'CD8+ TEMRA cells'
subset = subset(sc, condition %in% c('recovered_mild', 'recovered_severe') & integrated_annotations == cell.type & !!sym(dextramer) == 'YES')
Idents(subset) = 'condition'

markers = FindMarkers(subset, ident.1 = 'recovered_severe', ident.2 = 'recovered_mild', logfc.threshold = 0)
markers = markers[order(markers$avg_log2FC, decreasing = TRUE),]


# Write to file
markers.out = markers[markers$p_val_adj < 0.05,]
markers.out$gene = rownames(markers.out)
cell.type.name = gsub(' ', '_', cell.type)
file.name = paste0(cell.type.name, '_binding_recovered_severe_vs_recovered_mild')

write.table(markers.out[,c(6, 1:5)], 
			file = paste0(outdir, dextramer, '/', file.name, '_dge.txt'), 
			sep = '\t',
			row.names = FALSE, 
			quote = FALSE)


# GSEA
# GO
run_gsea(object = sc, markers = markers, category = 'C5', subcategory = 'BP',
		out.dir = paste0(outdir, dextramer), cell.type = cell.type, 
		file.prefix = file.name)

# PID
run_gsea(object = sc, markers = markers, category = 'C2', subcategory = 'PID',
		out.dir = paste0(outdir, dextramer), cell.type = cell.type, 
		file.prefix = file.name)

# Immunological signature
run_gsea(object = sc, markers = markers, category = 'C7',
		out.dir = paste0(outdir, dextramer), cell.type = cell.type,
		file.prefix = file.name)	









