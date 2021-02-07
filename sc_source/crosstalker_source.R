# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de


#---- Generate CrossTalkeR report

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




#---- Differential frequency table (modified from CrossTalkeR)
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




#---- Differential gene cell ranking (modified from CrossTalkeR)

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




#---- Sankey plots of interactions from CrossTalkeR

plot_sankey = function(data, target){
  data = data[grep(target, data$allpair),]
  data$Expression = ifelse(data$MeanLR > 0, 'Up', 'Down')
  data$AbsMeanLR = abs(data$MeanLR)
  
  p = ggplot(data, aes(y = AbsMeanLR,
                       axis1 = Ligand.Cluster,
                       axis2 = Ligand,
                       axis3 = Receptor,
                       axis4 = Receptor.Cluster)) +
    geom_alluvium(aes(fill = Expression), width = 1/12) +
    geom_stratum(width = 1/12) +
    geom_label(stat = 'stratum', aes(label = after_stat(stratum))) +
    scale_x_discrete(limits = c('Ligand.Cluster', 'Ligand', 
                                'Receptor', 'Receptor.Cluster')) +
    cowplot::theme_cowplot() +
    theme(axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.title.y = element_blank(), 
          axis.text.y = element_blank()) +
    scale_fill_manual(values = c('Up' = viridis(2)[2],
                      'Down' = viridis(2)[1]))
  
  return(p)
}

