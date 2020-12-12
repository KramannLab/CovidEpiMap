# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de


# Cell type coloring scheme
# Colors inspired by https://sashat.me/2017/01/11/list-of-20-simple-distinct-colors/
cell.type.colors = c('CD8+ naive T cells' = '#f58231',
					           'CD8+ central memory T cells' = '#ffe119',
                     'CD8+ CD73+ regulatory T cells' = '#fabebe',
					           'CD8+ TEMRA cells' = '#000075',
					           'CD8+ NK-like TEMRA cells' = '#568198',
					           'CD8+ effector memory T cells 1' = '#70036a',
					           'CD8+ effector memory T cells 2' = '#911eb4',
					           'CD8+ cycling effector T cells' = '#f032e6',
					           'CD8+ NK-like early effector T cells' = '#ff8ff9',
					           'Atypical NKT cells' = '#aaffc3',
					           'CD8+ exhausted T cells' = '#e6194B',
					           'MAIT cells' = '#999999',
					           'Gamma Delta T cells' = '#000000')




#---- Genesorter
run_genesorter = function(sc, assay = 'RNA', slot = 'data', write.file = FALSE, out.dir = '.', file.name.prefix = NULL){
  
  # Get specificity score (specScore) and conditional probability of expression (condGeneProb)
  library(genesorteR, quietly = TRUE)
  sg = sortGenes(GetAssayData(sc, assay = assay, slot = slot), Idents(sc))
  return.list = c('specScore' = sg$specScore, 'condGeneProb' = sg$condGeneProb)

  # Write files
  if(write.file){
    if(!dir.exists(file.path(out.dir))) stop('out.dir does not exist')

    # specScore
    specScore = as.data.frame(sg$specScore)
    specScore$gene = rownames(specScore)
    specScore = specScore[,c(ncol(specScore), 1:ncol(specScore)-1)]
    write.table(specScore, 
      file = paste0(out.dir, '/', file.name.prefix, 'specScore.txt'), 
      sep = '\t', quote = FALSE, row.names = FALSE)

    # condGeneProb
    condGeneProb = as.data.frame(sg$condGeneProb)
    condGeneProb$gene = rownames(condGeneProb)
    condGeneProb = condGeneProb[,c(ncol(condGeneProb), 1:ncol(condGeneProb)-1)]
    write.table(condGeneProb, 
      file =  paste0(out.dir, '/', file.name.prefix, 'condGeneProb.txt'), 
      sep = '\t', quote = FALSE, row.names = FALSE)
  }

  return(return.list)
}




#---- Correlation heatmap with meta data
correlation_heatmap = function(object, conditionVector, assay = 'RNA', cellTypeColors = cell.type.colors){
  suppressPackageStartupMessages(library(ComplexHeatmap))
  suppressPackageStartupMessages(library(stringr))
  suppressPackageStartupMessages(library(RColorBrewer))
  suppressPackageStartupMessages(library(viridis))
  DefaultAssay(object) = assay

  if(assay == 'RNA'){
      df = Reduce('rbind', 
              AverageExpression(object,
                features = VariableFeatures(object, assay = assay), 
              assays = assay))
  }

  if(assay == 'ADT'){
    # Exclude dextramers in correlation
      df = Reduce('rbind', 
              AverageExpression(object, 
                features = rownames(object)[23:38], assays = assay))
  }

  df.plot = cor(df)

  # Get meta data
  condition = str_extract(colnames(df), paste(conditionVector, collapse = '|'))
  cell.types = gsub(paste0('[[:space:]]?', conditionVector, collapse = '|'), '', colnames(df))

  # Heatmap annotations
  mat.colors = colorRampPalette(rev(brewer.pal(n = 7, name = 'RdYlBu')))(100)
  condition.colors = viridis(length(conditionVector))
  names(condition.colors) = conditionVector
  ha = HeatmapAnnotation(Condition = condition,
                        Celltypes = cell.types,
                        col = list(Condition = condition.colors,
                                  Celltypes = cellTypeColors),
                        show_annotation_name = FALSE)

  # Heatmap
  p = Heatmap(matrix = df.plot,
             col = mat.colors,
             name = 'Pearson',
             show_row_names = FALSE,
             show_column_names = FALSE,
             clustering_distance_rows = 'pearson',
             clustering_distance_columns = 'pearson',
             clustering_method_columns = 'average',
             clustering_method_rows = 'average',
             cluster_rows = cluster_within_group(df.plot, cell.types),
             cluster_columns = cluster_within_group(df.plot, cell.types),
             top_annotation = ha)
  return(p)
  rm(df)
  rm(df.plot)
}




#---- Differential gene expression function for integrated data
get_markers = function(sc, condition, control, cell.type){
  # Take care if no cells
  markers = data.frame()
  condition.subset = tryCatch(subset(sc, celltype.condition == paste(cell.type, condition)), 
    error = function(e) {data.frame()})
  control.subset = tryCatch(subset(sc, celltype.condition == paste(cell.type, control)), 
    error = function(e) {data.frame()})

  if (ncol(condition.subset) > 20 & ncol(control.subset) > 20){
    markers = FindMarkers(sc, ident.1 = paste(cell.type, condition), 
      ident.2 = paste(cell.type, control), min.pct = 0.25)
    markers = markers[markers$p_val_adj < 0.05,]
  
    markers$gene = rownames(markers)
    cols = colnames(markers)
    cell.type.name = gsub('/', '_',cell.type)
    cell.type.name = gsub(' ', '_',cell.type.name)

    if (nrow(markers) > 0){
      write.table(markers[,c('gene',cols[-6])], 
        file = paste0('diff_genes/integrated.diff.genes.', cell.type.name,
          '.', condition, '.vs.', control, '.txt'),
        quote = FALSE, sep = '\t', row.names = FALSE)
    }
  }
  return(markers)
}




#---- Heatmap function for DE genes of integrated data
top_heatmap = function(sc, markers, cell.type, disease, control, n){
  if (nrow(markers) > 0){
    top.genes = markers %>%
      top_n(n = n, wt = avg_log2FC)

    bot.genes = markers %>%
      top_n(n = -n, wt = avg_log2FC)

    genes = unique(c(top.genes$gene, bot.genes$gene))

    cell.type.name = gsub('/', '_',cell.type)
    cell.type.name = gsub(' ', '_',cell.type.name)
    pdf(file = paste0('diff_genes/heatmaps/integrated_diff_genes_', cell.type.name,
        '_', disease, '_vs_', control, '.pdf'))
    print(DoHeatmap(subset(sc, integrated_annotations == cell.type & condition %in% c(disease, control)), 
      features = genes, group.by = 'condition',
      size = 3, raster = FALSE, angle = 30) + labs(title = cell.type))
    dev.off()
  }
}


