# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de


#---- Per-sample clustering and integration of GEX and ADT libraries

# Source 
# https://satijalab.org/seurat/v4.0/weighted_nearest_neighbor_analysis.html

#remotes::install_github("satijalab/seurat", ref = "release/4.0.0")
library(Seurat)
library(SeuratData)
library(genesorteR)
library(stringr)
library(ggplot2)
library(dplyr)
source('sc_source/sc_source.R')
setwd('~/sciebo/CovidEpiMap/')
patients = c('1', '2', '3', '5', '6', '7', '11', '12', 
	'13', '15', '18', '19', '29', '31', '32')


#---- Multimodal integration per sample individually

for (patient in patients){
	# Define sample
	print(paste0('** Processing patient ', patient, ' **'))
	sample = paste0(patient, '_GEX_SURF')
	outdir = paste0('per_patient/', sample, '/')
	dir.create(outdir)



	# Read both Gene Expression and Antibody Capture (ADT) assays
	sc.data = Read10X(data.dir = paste0('data/', sample, '/filtered_feature_bc_matrix'))
	sc = CreateSeuratObject(counts = sc.data$'Gene Expression', project = sample)
	sc[['ADT']] = CreateAssayObject(counts = sc.data$'Antibody Capture')



	# QC plots before filtering
	DefaultAssay(sc) = 'RNA'
	sc[['percent.mt']] = PercentageFeatureSet(sc, pattern = '^MT-')
	sc[['percent.ribo']] = PercentageFeatureSet(sc, pattern ="^RP[SL]")

	features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt', 
			'percent.ribo', 'nFeature_ADT', 'nCount_ADT')

	pdf(file = paste0(outdir, sample, '_QC.pdf'), height = 10, width = 14)
	print(VlnPlot(sc, ncol = 3, pt.size = 0,
		features = features))
	dev.off()



	# Pre-processing and dimensional reduction of the assays independently
	# RNA
	sc = subset(sc, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 10)
	sc = NormalizeData(sc, normalization.method = 'LogNormalize', verbose = FALSE)
	sc = FindVariableFeatures(sc, nfeatures = 2000, verbose = FALSE)
	sc = ScaleData(sc, verbose = FALSE, features = rownames(sc))
	sc = RunPCA(sc, verbose = FALSE)



	# ADT
	# Do not use dextramer counts for dimensional reduction (first 23 features)
	DefaultAssay(sc) = 'ADT'
	VariableFeatures(sc) = rownames(sc[['ADT']])
	sc = NormalizeData(sc, normalization.method = 'CLR', margin = 2, verbose = FALSE)
	sc = ScaleData(sc, verbose = FALSE, features = rownames(sc))
	# Add another PCA name so not to overwrite RNA's PCA slot
	sc = RunPCA(sc, verbose = FALSE, reduction.name = 'apca', reduction.key = 'APC_',
		features = rownames(sc)[23:38])
	apca_dims = ncol(sc@reductions$apca@cell.embeddings)



	# Construct weighted nearest neighbors graph based on a weighted combination of RNA and ADT
	sc = FindMultiModalNeighbors(sc, reduction.list = list('pca', 'apca'),
		dims.list = list(1:30, 1:apca_dims), modality.weight.name = 'RNA.weight')



	# Visualize and cluster data
	# UMAP based on the weighted combination of RNA and ADT data
	sc = RunUMAP(sc, nn.name = 'weighted.nn', reduction.name = 'wnn.umap', reduction.key = 'wnnUMAP_')
	sc = FindClusters(sc, graph.name = 'wsnn', algorithm = 3, resolution = 2, verbose = FALSE)


	# Clustering based on GEX and ADT separately
	sc = RunUMAP(sc, reduction = 'pca', dims = 1:30, assay = 'RNA',
		reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
	sc = RunUMAP(sc, reduction = 'apca', dims = 1:apca_dims, assay = 'ADT',
		reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')


	pdf(file = paste0(outdir, sample, '_clustering.pdf'))
	print(DimPlot(sc, reduction = 'wnn.umap', label = TRUE))
	print(DimPlot(sc, reduction = 'rna.umap', label = TRUE, repel = TRUE) + NoLegend())
	print(DimPlot(sc, reduction = 'adt.umap', label = TRUE, repel = TRUE) + NoLegend())
	dev.off()



	# QC plots after filtering and clustering
	pdf(file = paste0(outdir, sample, '_QC_afterfiltering.pdf'), height = 10, width = 14)
	print(VlnPlot(sc, ncol = 3, pt.size = 0,
		features = features))
	print(FeaturePlot(sc, reduction = 'wnn.umap', ncol = 3, 
		features = features))
	dev.off()



	# Visualize counts from ADT library (dextramer + CITE seq)
	DefaultAssay(sc) = 'ADT'
	markers = rownames(sc[['ADT']])
	pdf(file = paste0(outdir, sample, '_ADT_Markers.pdf'), width = 8)
	print(DoHeatmap(sc, features = markers, assay = 'ADT', raster = FALSE, size = 3))
	for (marker in markers){
		print(FeaturePlot(sc, feature = marker, reduction = 'wnn.umap'))
		print(VlnPlot(sc, feature = marker))
	}
	dev.off()



	# Get gene specificity score and conditionial probability of expression per cluster
	# Use for cluster annotations
	sg = run_genesorter(sc, write.file = TRUE, file.name.prefix = sample)



	# Dotplot of top 10 genes per cluster based on specificity score
	sg = sg$specScore
	genes = list()
	n = 10
	for (i in 1:ncol(sg)){
		genes[[i]] = rownames(head(sg[order(sg[,i], decreasing = TRUE),], n))
	}
	top.genes = unique(unlist(genes))

	height = length(top.genes) * 0.2
	width = length(levels(Idents(sc))) * 0.5
	pdf(file = paste0(outdir, sample, '_GeneSorterTop10.pdf'),
		height = height, width = width)
	print(DotPlot(sc, features = rev(top.genes), assay = 'RNA') + coord_flip())
	dev.off()



	# Top 5 differentially expressed genes per cluster
	sc.markers = FindAllMarkers(sc, assay = 'RNA', only.pos = TRUE, min.pct = 0.25)
	sc.markers = sc.markers[sc.markers$p_val_adj < 0.05,]
	cols = colnames(sc.markers)
	write.table(sc.markers[,c('gene',cols[-7])], 
		file = paste0(outdir, sample, '.diff.genes.txt'),
		quote = FALSE, sep = '\t', row.names = FALSE)

	top.genes = sc.markers %>%
		group_by(cluster) %>%
		top_n(n = 5, wt = avg_log2FC)

	width = length(levels(Idents(sc))) * 0.8
	height = length(unique(top.genes$gene)) * 0.2
	pdf(file = paste0(outdir, sample, '_Top5Markers.pdf'), 
		width = width, height = height)
	print(DoHeatmap(sc, features = unique(top.genes$gene), size = 3,
		raster = FALSE, angle = 90) + NoLegend())
	dev.off()



	# Save data
	DefaultAssay(sc) = 'RNA'
	saveRDS(sc, file = paste0(outdir, sample, '.rds'))

}



#---- Add TCR clonotypes to Seurat objects

# Source
# https://www.biostars.org/p/383217/

# Run addVDJ.py script prior to the following

for (patient in patients){

	# Define sample
	sample = paste0(patient, '_GEX_SURF')
	indir = paste0('per_patient/', sample, '/')
	sc = readRDS(file = paste0(indir, sample, '.rds'))



	# Add TCR clonotypes
	vdj_type = 'TCR'
	
	# Read VDJ info
	sample_vdj = paste0(patient, '_', vdj_type)
	vdj = read.csv(file = paste0('data/', sample_vdj, '/filtered_contig_annotations.csv'))

	# Subset so only first line of each barcode is kept (they have the same clonotype id anyway)
	vdj = vdj[!duplicated(vdj$barcode),]
	vdj = vdj[,c('barcode', 'raw_clonotype_id')]
	names(vdj)[names(vdj) == 'raw_clonotype_id'] = 'clonotype_id'

	# Get clonotype info
	clono = read.csv(file = paste0('data/', sample_vdj, '/clonotypes.csv'))

	# Add amino acid sequence of CDR3s to data frame
	vdj = merge(vdj, clono[,c('clonotype_id', 'cdr3s_aa')])

	# Re-order to get barcodes as first column again
	vdj = vdj[,c(2,1,3)]
	rownames(vdj) = vdj[,1]
	vdj[,1] = NULL
	colnames(vdj) = c(paste0(vdj_type, '_clonotype_id'), paste0(vdj_type, '_cdr3s_aa'))
	cells = rownames(vdj)



	# Add VDJC genes
	vdj.genes = read.csv(file = paste0(sample_vdj, '/clonotype.vdj.genes.csv'))
	colnames(vdj.genes) = sub('Clonotype', paste0(vdj_type, '_clonotype_id'), colnames(vdj.genes))
	vdj = merge(vdj, vdj.genes, paste0(vdj_type, '_clonotype_id'))
	rownames(vdj) = cells



	# Add clonotype information to Seurat object
	sc = AddMetaData(object = sc, metadata = vdj)
	saveRDS(sc, file = paste0(indir, sample, '.rds'))
}


