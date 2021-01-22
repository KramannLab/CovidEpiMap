# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de


#---- Cell type coloring scheme

# Source 
# https://sashat.me/2017/01/11/list-of-20-simple-distinct-colors/
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

get_markers = function(sc, condition, control, cell.type, only.sig = TRUE, adj.pval.cutoff = 0.05, logfc.threshold = 0.25, out.dir = 'diff_genes'){
  # Take care if no cells
  markers = data.frame()
  condition.subset = tryCatch(subset(sc, celltype.condition == paste(cell.type, condition)), 
    error = function(e) {data.frame()})
  control.subset = tryCatch(subset(sc, celltype.condition == paste(cell.type, control)), 
    error = function(e) {data.frame()})

  if (ncol(condition.subset) > 20 & ncol(control.subset) > 20){
    markers = FindMarkers(sc, ident.1 = paste(cell.type, condition), 
      ident.2 = paste(cell.type, control), min.pct = 0.25, logfc.threshold = logfc.threshold)

    if (only.sig){
      markers = markers[markers$p_val_adj < adj.pval.cutoff,]
    }
    
  
    markers$gene = rownames(markers)
    cols = colnames(markers)
    cell.type.name = gsub('/', '_',cell.type)
    cell.type.name = gsub(' ', '_',cell.type.name)

    if (nrow(markers) > 0){
      write.table(markers[,c('gene',cols[-6])], 
        file = paste0(out.dir, '/integrated.diff.genes.', cell.type.name,
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




#---- Cell-type abundance testing

# Source
# https://www.nxn.se/valent/2020/11/28/s9jjv32ogiplagwx8xrkjk532p7k28
# Uses integrated annotations!
get_abundance = function(df, condition.vector){
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(tidyverse))
  suppressPackageStartupMessages(library(ggrepel))
  suppressPackageStartupMessages(library(emmeans))
  if(length(condition.vector) > 2) stop('Only two conditions allowed')

  # Prepare linear model
  data.df = df %>% filter(condition %in% condition.vector)
  formula = cbind(count, other) ~ integrated_annotations * condition
  model1 = glm(formula = formula, family = 'binomial', data = data.df)


  # Log odds ratios
  emm1 = emmeans(model1, specs = pairwise ~ condition | integrated_annotations)
  c_results = emm1$contrasts %>%
  summary(infer = TRUE, type = 'response') %>%
    rbind() %>%
    as.data.frame()


  p1 = (
    ggplot(data = c_results, aes(x = odds.ratio, y = -log10(p.value))) +
      geom_point() + 
      geom_text_repel(data = c_results %>% filter(p.value < 1e-50), 
            aes(label = integrated_annotations)) +
      scale_x_log10() +
      theme_minimal() +
        labs(x = paste0(condition.vector[2], ' / ', condition.vector[1]),
          title = 'Cell type proportion change')
    )


  # Relate log odds ratios with the average abundace of each cell type
  emm2 = emmeans(model1, specs = ~ integrated_annotations)
  mean_probs = emm2 %>%
  summary(type = 'response') %>%
  select(integrated_annotations, prob)

  m_results = c_results %>% left_join(mean_probs)


  p2 = (
      ggplot(data = m_results, aes(x = prob, y = odds.ratio, color = p.value < 0.001)) +
        geom_point() +
        geom_text_repel(data = m_results %>% filter(abs(log(odds.ratio)) > log(2)), 
                aes(label = integrated_annotations), 
                color = 'black') +
        scale_x_log10() +
        scale_y_log10() +
        theme_minimal() +
        labs(y = paste0(condition.vector[2], ' / ', condition.vector[1], ' (odds ratio)'), 
          title = 'Cell type proportion change', 
          x = 'Average abundance (probability)')
  )

  return(list(p1 = p1, p2 = p2))

  }




#---- TF activity inference with DoRothEA (regulon has to be built prior to this)

run_dorothea = function(case, control, diff.indir, out.dir, dorothea_regulon_human, regulon){
  suppressPackageStartupMessages(library(viper))

  TF_activities_df = data.frame(tf = unique(dorothea_regulon_human$tf), 
              row.names = unique(dorothea_regulon_human$tf))

  # Get cell types
  pattern = paste0('.', case, '.vs.', control, '.txt')
  case.files = list.files(diff.indir, pattern = pattern)
  cell.types = sub('integrated.diff.genes.', '', case.files)
  cell.types = sub(pattern, '', cell.types)


  # Run DoRothEA on case/control contrasts per cell type
  for (cell.type in cell.types){
    file.name = paste0(diff.indir, 'integrated.diff.genes.', 
          gsub(' ', '_', cell.type), '.', 
          case, '.vs.', control, '.txt')

    diff.genes = read.table(file.name, header = TRUE)

    # Estimate z-score values for the gene expression signature (GES)
    myStatistics = matrix(diff.genes$avg_log2FC, dimnames = list(diff.genes$gene, 'avg_log2FC'))
    myPvalue = matrix(diff.genes$p_val, dimnames = list(diff.genes$gene, 'p_val'))
    mySignature = (qnorm(myPvalue/2, lower.tail = FALSE) * sign(myStatistics))[, 1]
    mySignature = mySignature[order(mySignature, decreasing = TRUE)]

    # Estimate TF activity
        mrs = msviper(ges = mySignature, regulon = regulon, minsize = 4, 
              ges.filter = FALSE, verbose = FALSE)
        TF_activities = data.frame(Regulon = names(mrs$es$nes),
                           Size = mrs$es$size[ names(mrs$es$nes) ], 
                           NES = mrs$es$nes, 
                           p.value = mrs$es$p.value, 
                           FDR = p.adjust(mrs$es$p.value, method = 'fdr'))
        TF_activities = TF_activities[order(TF_activities$FDR),]

    # Write to file per cell type
    write.table(TF_activities, 
      file = paste0(out.dir, case, '_vs_', control, '_tf_activity_', 
                    gsub(' ', '_', cell.type), '.txt'),
      sep = '\t',
      row.names = FALSE,
      quote = FALSE)

    TF_activities_df[[cell.type]] = TF_activities[rownames(TF_activities_df),'NES']
  }
  return(TF_activities_df)
}




#---- Plot results from run_dorothea

plot_dorothea = function(df, case, control){
  suppressPackageStartupMessages(library(pheatmap))

  # Order heatmap by TF activity
  df$tf = NULL
  df.plot = df[complete.cases(df),]
  tf_cluster = hclust(dist(t(df.plot)))
  tfs_ordered = tf_cluster$labels[tf_cluster$order]
  celltype_cluster = hclust(dist(df.plot))
  celltype_ordered = celltype_cluster$labels[celltype_cluster$order]

  paletteLength = 100
  myColor = colorRampPalette(c('Darkblue', 'white', 'red'))(paletteLength)
  viperBreaks = c(seq(min(df.plot), 0, 
                      length.out = ceiling(paletteLength/2) + 1),
                  seq(max(df.plot)/paletteLength, 
                      max(df.plot), 
                      length.out = floor(paletteLength/2)))

  p = pheatmap(t(df.plot)[tfs_ordered,celltype_ordered],
        fontsize_col = 4, 
        fontsize_row = 8,
        color = myColor, 
        breaks = viperBreaks, 
        main = paste0(case, ' vs ', control),
        angle_col = 90,
        border_color = NA)

  return(p)
}




#---- Plot GO terms from GSEA

plot_go = function(gsea.res, gsea.res.order, plot.title = NULL, n = 20){
  suppressPackageStartupMessages(library(ggplot2))

  # Format GOs
  plot.table = head(gsea.res[gsea.res.order,], n = n)
  plot.table$pathway = sub('GO_', '', plot.table$pathway)
  plot.table$pathway = gsub('_', ' ', plot.table$pathway)

  p = ggplot(plot.table,
    aes(x = NES, y = pathway)) +
    geom_point(aes(size = NES, color = padj)) +
    theme_bw(base_size = 8) +
    ylab(NULL) +
    ggtitle(plot.title) +
    scale_colour_gradient2(low = 'red', 
                mid = 'lightgrey', 
                high = 'blue', 
                midpoint = 0.05, 
                limits = c(0,0.1), 
                oob = scales::squish)

  return(p)
}




#---- GSEA function

run_gsea = function(bg.genes, stats, category, plot.title = NULL, subcategory = NULL, out.dir = '.', file.prefix, n = 30){
  suppressPackageStartupMessages(library(msigdbr))
  suppressPackageStartupMessages(library(fgsea))
  suppressPackageStartupMessages(library(dplyr))

  # Fetch geneset
  geneSets = msigdbr(species = 'Homo sapiens', category = category, subcategory = subcategory)
  geneSets = geneSets[geneSets$human_gene_symbol %in% bg.genes,]
  m_list = geneSets %>% split(x = .$human_gene_symbol, f = .$gs_name)

  # Run GSEA
  gsea = fgsea(pathways = m_list, stats = stats, minSize = 10, eps = 0.0)
  order = order(gsea$padj, decreasing = FALSE)

  # Plot
  file.name = paste0(out.dir, '/', file.prefix, '_gsea_', paste0(c(category, subcategory), collapse = '_'))

  pdf(file = paste0(file.name, '.pdf'), width = 10, height = 9)
  print(plot_go(gsea.res = gsea,
      gsea.res.order = order, 
      n = n, 
      plot.title = plot.title))
  dev.off()

  # Write to file
  write.table(gsea[order, -8], 
      file = paste0(file.name, '.txt'), 
      sep = '\t', 
      row.names = FALSE, 
      quote = FALSE)
}




#---- Nice GO plot

plot_go_nice = function(gse, terms, title = NULL){
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(viridis))
  
  # Format terms 
  terms.tmp = stringr::str_to_upper(terms)
  terms.tmp = gsub(' ', '_', terms.tmp)
  
  # Subset table
  gse = gse[gse$pathway %in% terms.tmp,]
  rownames(gse) = gse$pathway
  gse = gse[terms.tmp,]
  gse$pathway = terms
  gse$size = abs(gse$NES)
  
  # Plot
  p = gse %>%
    ggplot(aes(x = NES, y = pathway, color = padj, size = size)) +
    geom_point() +
    cowplot::theme_cowplot() + 
    theme(axis.ticks = element_blank(),
          plot.title = element_text(hjust = 0, size = 4),
          axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size = 8),
          axis.title.x = element_text(size = 10)) +
    scale_colour_gradient2(low = viridis(2)[2], 
                           mid = 'lightgrey', 
                           high = viridis(2)[1], 
                           midpoint = 0.05, 
                           limits = c(0,0.1), 
                           oob = scales::squish,
                           name = 'pAdj') +
    labs(x = 'NES',
         y = '',
         size = '',
         title = title) +
    guides(size = 'none') +
    scale_size(range = c(1,5),
               limits = c(1,3)) +
    coord_fixed(xlim = c(-3.5, 3.5))
  return(p)
}




#---- Plot smoothed expression of genes from tradeSeq

plot_top_genes = function(gene.set, model, n, out.dir, file.name){
  suppressPackageStartupMessages(library(tradeSeq))
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(SingleCellExperiment))

  pdf(file = paste0(out.dir, '/', file.name, '.pdf'))
  for (gene in gene.set[1:n]){
    print(plotSmoothers(model, assays(model)$counts,
      gene = gene,
      alpha = 1,
      border = TRUE) +
      ggtitle(gene))
  }
  dev.off()
}




#----Combined marker violin plot

get_violin = function(object, features.use){
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(Seurat))
  suppressPackageStartupMessages(library(cowplot))

    p = VlnPlot(object = object, 
                 features = features.use, 
                 pt.size = 0, cols = cell.type.colors,
                 combine = FALSE)
    
    ps = lapply(p, function(x) x + coord_flip() + NoLegend() +
                     theme_bw() +
                     theme(plot.title = element_text(angle = 90, size = 8, vjust = 0.5), 
                           legend.position = 'none',
                           axis.title.x = element_blank(),
                           axis.text.x = element_blank(),
                           axis.ticks.x = element_blank(),
                           axis.text.y = element_blank(),
                           axis.title.y = element_blank(),
                           axis.ticks.y = element_blank(),
                           plot.margin = unit(c(0, 0, 0, 0), 'cm'),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           panel.background = element_rect(colour = 'black', size = 1)))
    
    p = plot_grid(plotlist = ps, nrow = 1, align = 'h')
    
    return(p)
}




#---- Compute p-values and make violin plots for Progeny pathway inference results

compute_stats = function(df, celltype, pathways, conditions, plot = FALSE){
  suppressPackageStartupMessages(library(tidyverse))
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(viridis))
  
  p = NULL
  df_sub = df %>% filter(condition %in% conditions,
                         Pathway %in% pathways,
                         cell.type %in% celltype)
  
  # p.adjust for number of cells * number of pathways tested
  num_test = nrow(df_sub) * length(pathways)
  
  stats = df_sub %>% group_by(Pathway) %>% 
    nest() %>%
    mutate(wilcox = map(data, function(df){
      stest = wilcox.test(Activity ~ condition, data = df, alternative = 'two.sided')
      broom::tidy(stest) %>%
        mutate(corr_pvalue = p.adjust(p.value, method = 'BH', n = num_test))})) %>%
    select(-data) %>%
    unnest(wilcox) %>%
    ungroup() %>%
    arrange(corr_pvalue) %>% 
    as.data.frame
  
  if(plot){
    p = ggplot(df_sub, aes(x = condition, y = Activity, fill = condition)) +
      geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), colour = 'lightgrey') +
      theme_classic() + 
      scale_fill_viridis(discrete = TRUE, option = 'viridis') +
      theme(strip.background = element_rect(fill = 'lightgrey'),
            axis.ticks = element_blank()) +
      xlab('') +
      ylab('Pathway activity') +
      facet_wrap( ~ Pathway)
  }
  return(list('stats' = stats, 'p' = p))
}




#---- Nice dot plot of differential expression results from FindMarkers

# Source
# https://davemcg.github.io/post/lets-plot-scrna-dotplots/
plot_dge_nice = function(dge.table, genes){
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(viridis))

  # Get genes and stats to plot
  gene.list = list()
  for (i in 1:length(dge.table)){
    table = dge.table[[i]]
    
    col = table %>% 
      filter(gene %in% genes) %>% 
      select(gene, avg_log2FC, p_val_adj) %>% 
      mutate(cluster = names(dge.table[i]),
             log10_padj = -log10(p_val_adj))
    
    gene.list[[i]] = col
  }
  gene.table = do.call(rbind, gene.list)
  gene.table$log10_padj_edit = gene.table$log10_padj
  gene.table[gene.table$log10_padj > 15,'log10_padj_edit'] = 15
  gene.table$gene = factor(gene.table$gene, levels = rev(genes))
  
  # Plot
  p = gene.table %>%
    ggplot(aes(x = cluster, y = gene, color = avg_log2FC, size = log10_padj_edit)) +
    geom_point() +
    scale_colour_gradient2(low = viridis(2)[1], 
                           mid = '#F0F0F0', 
                           high = viridis(2)[2], 
                           midpoint = 0, 
                           limits = c(-2,2), 
                           oob = scales::squish) +
    cowplot::theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
          axis.text.y = element_text(size = 8),
          axis.ticks = element_blank()) +
    labs(x = '',
         y = '',
         size = bquote(~-log[10]~'(pAdj)'),
         color = bquote(~log[2]~'(FC)'))
  
  return(p)
}




#---- Discrete colour gradient function 

# Source 
# https://stackoverflow.com/a/62556763
discrete_gradient_pal = function(colours, bins = 5) {
  ramp = scales::colour_ramp(colours)
  function(x) {
    if (length(x) == 0) return(character())
    i = floor(x * bins)
    i =  ifelse(i > bins-1, bins-1, i)
    ramp(i/(bins-1))
  }
}




#---- # Get upper triangle of matrix

get_upper_triangle = function(cormat, diag = FALSE){
  cormat[lower.tri(cormat, !diag)] =  NA
  return(cormat)
}




#---- Plot upper triangular matrix as a tile plot

upper_trig_tile_plot = function(matrix, text_size = 2){
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(cowplot))
  suppressPackageStartupMessages(library(viridis))
  colours = viridis(7)
  
  p = ggplot(coef_matrix, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() + 
    geom_text(aes(label = round(value, digits = 2)), size = text_size) +
    scale_fill_gradient2(high = colours[7], 
                         mid = colours[4], 
                         midpoint = ((range(na.omit(coef_matrix$value)))/2)[2], 
                         low = colours[1], 
                         na.value = 'white') +
    theme_cowplot() + 
    theme(text = element_text(size = 8),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_text(size = 8), 
          axis.text.x = element_text(size = 8, angle = 90, vjust = 1, hjust = 1))
  
  return(p)
}

