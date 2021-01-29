# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de


#---- Visualize CellPhoneDB results with CrossTalkeR

library(clusterProfiler)
library(org.Hs.eg.db)
library(patchwork)
library(dplyr)
library(viridis)
source('sc_source/sc_source.R')
source.folder = '~/sciebo/CovidEpiMap/CrossTalkeR/CrossTalkeR/'
devtools::load_all(source.folder)
wkdir = '~/sciebo/CovidEpiMap/CrossTalkeR/'
indir = 'significant_means/filtered_corrected/'
setwd(wkdir)


crosstalker_report = function(Control, Case, Genes, Indir, Outdir){
	# Define path for control and case (CellPhoneDB has been run separately on these)
	# Files have already been filtered and corrected with format.py from CrossTalkeR package
	paths = c('CTR' = paste0(Indir, Control, '_filtered_corrected.csv'),
			'EXP' = paste0(Indir, Case, '_filtered_corrected.csv'))

	# Generate report
	outdir = paste0(Outdir, '/', Control, '_vs_', Case, '/')

	generate_report(lrpaths = paths,
					out_path = outdir,
					threshold = 0,
					out_file = paste0(Control, '_vs_', Case, '.html'),
					report = TRUE)
}


# Healthy vs active mild
control = 'healthy'
case = 'active_mild'
crosstalker_report(Control = control, Case = case, Indir = indir, Outdir = wkdir)


# Healthy vs active severe
control = 'healthy'
case = 'active_severe'
crosstalker_report(Control = control, Case = case, Indir = indir, Outdir = wkdir)


# Active mild vs active severe
control = 'active_mild'
case = 'active_severe'
data = crosstalker_report(Control = control, Case = case, Indir = indir, Outdir = wkdir)


# Recovered mild vs recovered severe
control = 'recovered_mild'
case = 'recovered_severe'
crosstalker_report(Control = control, Case = case, Indir = indir, Outdir = wkdir)



#---- Save selected plots from CrossTalkeR (active severe vs active mild)

outdir = paste0(wkdir, 'active_mild_vs_active_severe/')


# Single cell-cell interaction plot (active mild)
pdf(file = paste0(outdir, 'single_cci_plot_active_mild.pdf'), width = 10, height = 10)
plot_cci(graph = data@graphs$CTR,
         colors = cell.type.colors,
         plt_name = '',
         coords = data@coords[V(data@graphs$CTR)$name,],
         emax = NULL,
         leg = TRUE,
         low = 0,
         high = 0,
         ignore_alpha = FALSE,
         log = FALSE,
         efactor = 8,
         vfactor = 12)
dev.off()

# Single cell-cell interaction plot (active severe)
pdf(file = paste0(outdir, 'single_cci_plot_active_severe.pdf'), width = 10, height = 10)
plot_cci(graph = data@graphs$EXP,
         colors = cell.type.colors,
         plt_name = '',
         coords = data@coords[V(data@graphs$EXP)$name,],
         emax = NULL,
         leg = TRUE,
         low = 0,
         high = 0,
         ignore_alpha = FALSE,
         log = FALSE,
         efactor = 8,
         vfactor = 12)
dev.off()

# Differential cell-cell interaction plot
pdf(file = paste0(outdir, 'differential_cci_plot.pdf'), width = 10, height = 10)
plot_cci(graph = data@graphs$EXP_x_CTR,
         colors = cell.type.colors,
         plt_name = '',
         coords = data@coords[V(data@graphs$EXP_x_CTR)$name,],
         emax = NULL,
         leg = TRUE,
         low = 0,
         high = 0,
         ignore_alpha = FALSE,
         log = FALSE,
         efactor = 8,
         vfactor = 12)
dev.off()



# Plot differential cell-cell ligand/receptor frequency table
all_data = readRDS(file = paste0(outdir, 'LR_data_step2.Rds'))

# Plotting script from CrossTalkeR (modified)
plot_table = function(all_data = all_data){
    curr <- "EXP_x_CTR"
    curr_net <- all_data@graphs[[curr]]
    in_deg_up <- table(all_data@tables[[curr]]$Ligand.Cluster[all_data@tables[[curr]]$MeanLR > 0])
    in_up <- tibble::tibble(as.data.frame(in_deg_up))
    in_deg_down <- table(all_data@tables[[curr]]$Ligand.Cluster[all_data@tables[[curr]]$MeanLR < 0])
    in_down <- tibble::tibble(as.data.frame(in_deg_down))
    in_down$Freq <- 0-in_down$Freq
    in_all <- dplyr::bind_rows(in_up,in_down)
    in_all$Expression <- ifelse(in_all$Freq<0,'Downregulated','Upregulated')
    in_all$rank <- ifelse(in_all$Freq<0,0,1)
    out_deg_up <- table(all_data@tables[[curr]]$Receptor.Cluster[all_data@tables[[curr]]$MeanLR > 0])
    out_up <- tibble::tibble(as.data.frame(out_deg_up))
    out_deg_down <- table(all_data@tables[[curr]]$Receptor.Cluster[all_data@tables[[curr]]$MeanLR < 0])
    out_down <- tibble::tibble(as.data.frame(out_deg_down))
    out_down$Freq <- 0-out_down$Freq
    out_all <- dplyr::bind_rows(out_up,out_down)
    out_all$Expression <- ifelse(out_all$Freq<0,'Downregulated','Upregulated')
    out_all$rank <- ifelse(out_all$Freq<0,0,1)
    p1 <- ggplot2::ggplot(in_all,ggplot2::aes(x=Freq,y=Var1,fill=Expression))+
      ggplot2::geom_bar(stat = 'identity',position = "identity", show.legend = FALSE)+
      ggplot2::geom_text(ggplot2::aes(label=sub("-", "", Freq)),size=3.5)+
      ggplot2::scale_fill_manual(values=pals::coolwarm(2))+
      ggplot2::ggtitle('Ligands')+
      ggplot2::ylab('')+
      ggplot2::xlab('Number of interactions')+
      cowplot::theme_cowplot() +
      theme(axis.line  = element_blank(),
            axis.ticks = element_blank(),
            axis.text.x = element_blank())
    p2 <- ggplot2::ggplot(out_all,ggplot2::aes(x=Freq,y=Var1,fill=Expression))+
      ggplot2::geom_bar(stat = 'identity',position = "identity")+
      ggplot2::geom_text(ggplot2::aes(label=sub("-", "", Freq)),size=3.5)+
      ggplot2::scale_fill_manual(values=pals::coolwarm(2))+
      ggplot2::ggtitle('Receptors')+
      ggplot2::ylab('')+
      ggplot2::xlab('Number of interactions')+
      cowplot::theme_cowplot() +
      theme(axis.line  = element_blank(),
            axis.ticks = element_blank(),
            axis.text.x = element_blank())
    return(p1+p2)
}


pdf(file = paste0(outdir, 'differential_cci_frequency_table.pdf'), width = 13)
plot_table(all_data)
dev.off()



# Plot differential gene cell ranking (sending/receiving)
# Plotting script from CrossTalkeR (modified)
plot_rank = function(all_data){
  curr <- "EXP_x_CTR"
  curr_net <- all_data@graphs_ggi[[curr]]
  up_graph <- igraph::subgraph.edges(curr_net, E(curr_net)[E(curr_net)$MeanLR > 0])
  down_graph <- igraph::subgraph.edges(curr_net, E(curr_net)[E(curr_net)$MeanLR < 0])
  in_deg_up <- igraph::degree(up_graph, mode = 'in')
  in_deg_down <- igraph::degree(down_graph, mode = 'in')
  in_up <- tibble::tibble(genes = paste0(names(in_deg_up),'_up'), values=as.array(in_deg_up))
  in_down <- tibble::tibble(genes = paste0(names(in_deg_down),'_down'), values=as.array(in_deg_down))
  in_deg_data_up <-dplyr::top_n(in_up, 10, values)
  in_deg_data_down <- dplyr::top_n(in_down, 10, values)
  in_deg_data_down$values <- 0 -in_deg_data_down$values
  in_deg_data <- dplyr::bind_rows(in_deg_data_up,in_deg_data_down )
  in_deg_data$Expression <- ifelse(in_deg_data$values <0,'Downregulated','Upregulated')
  p1 <- ggplot2::ggplot(in_deg_data,ggplot2::aes(x=values,y=reorder(genes,values),fill=Expression))+
    ggplot2::geom_bar(stat = 'identity',position = "identity")+
    ggplot2::scale_fill_manual(values=pals::coolwarm(2))+
    ggplot2::geom_text(ggplot2::aes(label=sub("-", "",values)),size=3.5)+
    ggplot2::ggtitle('Receiving')+
    ggplot2::ylab('')+
    ggplot2::xlab('Number of interactions')+
    cowplot::theme_cowplot() +
    theme(axis.line  = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_blank())
  out_deg_up <- igraph::degree(up_graph, mode = 'out')
  out_deg_down <- igraph::degree(down_graph, mode = 'out')
  out_up <- tibble::tibble(genes = paste0(names(out_deg_up),'_up'), values=as.array(out_deg_up))
  out_down <- tibble::tibble(genes = paste0(names(out_deg_down),'_down'), values=as.array(out_deg_down))
  out_deg_data_up <-dplyr::top_n(out_up, 10, values)
  out_deg_data_down <- dplyr::top_n(out_down, 10, values)
  out_deg_data_down$values <- 0-out_deg_data_down$values
  out_deg_data <- dplyr::bind_rows(out_deg_data_up,out_deg_data_down )
  out_deg_data$Expression <- ifelse(out_deg_data$values <0,'Downregulated','Upregulated')
  p2 <- ggplot2::ggplot(out_deg_data,ggplot2::aes(x=values,y=reorder(genes,values),fill=Expression))+
    ggplot2::geom_bar(stat = 'identity',position = "identity", show.legend = FALSE)+
    ggplot2::scale_fill_manual(values=pals::coolwarm(2))+
    ggplot2::geom_text(ggplot2::aes(label=sub("-", "",values)),size=3.5)+
    ggplot2::ggtitle('Sending')+
    ggplot2::ylab('')+
    ggplot2::xlab('Number of interactions')+
    cowplot::theme_cowplot() +
    theme(axis.line  = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_blank())
  return(p2+p1)
}


pdf(file = paste0(outdir, 'differential_gene_cell_ranking.pdf'), width = 14)
plot_rank(all_data)
dev.off()



#---- KEGG pathway analysis

outdir = 'recovered_mild_vs_recovered_severe/'
data = readRDS(file = paste0(outdir, 'LR_data_final.Rds'))
cell.types = unique(c(unname(data@tables$EXP_x_CTR$Ligand.Cluster),
  unname(data@tables$EXP_x_CTR$Receptor.Cluster)))


pdf(file = paste0(outdir, 'KEGG_enrichment_significant.pdf'), height = 5, width = 7)
for (cell.type in cell.types){
  # Get interactions up/down in active mild vs active severe
  up = data@tables$EXP_x_CTR %>%
    filter(MeanLR > 0 & Ligand.Cluster %in%  cell.type)
  down = data@tables$EXP_x_CTR %>%
    filter(MeanLR < 0 & Ligand.Cluster %in%  cell.type)
  
  # Get unique up/down interactions
  up_exclu =  unique(up$Ligand)[match(unique(up$Ligand),unique(down$Ligand), nomatch = -1) == -1]
  down_exclu = unique(down$Ligand)[match(unique(down$Ligand),unique(up$Ligand), nomatch = -1) == -1]
  
  # KEGG enrichment analysis
  genes_up =  bitr(up_exclu, fromType = 'SYMBOL', toType = c('ENTREZID', 'ENSEMBL'), OrgDb = org.Hs.eg.db)
  genes_dw =  bitr(down_exclu, fromType = 'SYMBOL', toType = c('ENTREZID', 'ENSEMBL'), OrgDb = org.Hs.eg.db)
  enrich_up = enrichKEGG(genes_up$ENTREZID, organism = 'hsa')
  enrich_dw = enrichKEGG(genes_dw$ENTREZID, organism = 'hsa')
  
  # Plot significant results
  plot_data = enrich_up@result %>%
    mutate(log10.p.adjust = -log10(p.adjust)) %>%
    filter(p.adjust < 0.05) %>%
    top_n(n = 50, wt = log10.p.adjust)
  
  p_up = ggplot(plot_data, aes(x = log10.p.adjust, 
                          y = reorder(Description,log10.p.adjust))) +
    geom_bar(stat = 'identity', fill = viridis(2)[2]) +
    cowplot::theme_cowplot() +
    theme(axis.ticks = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          title = element_text(size = 12)) +
    labs(x = bquote(~-Log[10]~'(pAdj)'),
         title = paste0(cell.type, ': Up interactions'))
  
  plot_data = enrich_dw@result %>%
    mutate(log10.p.adjust = -log10(p.adjust)) %>%
    filter(p.adjust < 0.05) %>%
    top_n(n = 50, wt = log10.p.adjust)
  
  p_dw = ggplot(plot_data, aes(x = log10.p.adjust, 
                          y = reorder(Description,log10.p.adjust))) +
    geom_bar(stat = 'identity', fill = viridis(2)[1]) +
    cowplot::theme_cowplot() +
    theme(axis.ticks = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          title = element_text(size = 12)) +
    labs(x = bquote(~-Log[10]~'(pAdj)'),
         title = paste0(cell.type, ': Down interactions'))
  
  print(p_up)
  print(p_dw)
}
dev.off()

