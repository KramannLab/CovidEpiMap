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



