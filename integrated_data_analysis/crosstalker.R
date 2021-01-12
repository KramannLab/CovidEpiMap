# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de


#---- Visualize CellPhoneDB results with CrossTalkeR

source('sc_source/sc_source.R')
source.folder = '~/sciebo/CovidEpiMap/CrossTalkeR/CrossTalkeR/'
devtools::load_all(source.folder)
wkdir = '~/sciebo/CovidEpiMap/CrossTalkeR/'
indir = 'significant_means/filtered_corrected/'
setwd(wkdir)


crosstalker_report = function(Control, Case, Genes, Indir, Outdir){
	# Define path for control and case (CellPhoneDB has been run separately on these)
	# Files have already been filtered and corrected with organize.py from CrossTalkeR package
	paths = c('CTR' = paste0(Indir, Control, '_filtered_corrected.csv'),
			'EXP' = paste0(Indir, Case, '_filtered_corrected.csv'))

	# Generate report
	outdir = paste0(Outdir, '/', Control, '_vs_', Case, '/')

	generate_report(lrpaths = paths,
					genes = Genes,
					out_path = outdir,
					threshold = 0,
					out_file = paste0(Control, '_vs_', Case, '.html'),
					report = TRUE)
}


# Healthy vs active mild
control = 'healthy'
case = 'active_mild'
genes = c('MIF', 'KLRK1')
crosstalker_report(Control = control, Case = case, Genes = genes, Indir = indir, Outdir = wkdir)


# Healthy vs active severe
control = 'healthy'
case = 'active_severe'
genes = c('MIF', 'KLRK1')
crosstalker_report(Control = control, Case = case, Genes = genes, Indir = indir, Outdir = wkdir)


# Healthy vs recovered mild
control = 'healthy'
case = 'recovered_mild'
genes = c('MIF', 'KLRK1')
crosstalker_report(Control = control, Case = case, Genes = genes, Indir = indir, Outdir = wkdir)


# Healthy vs recovered severe
control = 'healthy'
case = 'recovered_severe'
genes = c('MIF', 'KLRK1')
crosstalker_report(Control = control, Case = case, Genes = genes, Indir = indir, Outdir = wkdir)


# Active mild vs active severe
control = 'active_mild'
case = 'active_severe'
genes = c('MIF', 'KLRK1')
data = crosstalker_report(Control = control, Case = case, Genes = genes, Indir = indir, Outdir = wkdir)


# Recovered mild vs recovered severe
control = 'recovered_mild'
case = 'recovered_severe'
genes = c('MIF', 'KLRK1')
crosstalker_report(Control = control, Case = case, Genes = genes, Indir = indir, Outdir = wkdir)



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
      ggplot2::geom_text(ggplot2::aes(label=Freq),size=3.5)+
      ggplot2::scale_fill_manual(values=pals::coolwarm(2))+
      ggplot2::ggtitle('Ligands')+
      ggplot2::ylab('')+
      ggplot2::xlab('Number of interactions')+
      ggplot2::theme_minimal()
    p2 <- ggplot2::ggplot(out_all,ggplot2::aes(x=Freq,y=Var1,fill=Expression))+
      ggplot2::geom_bar(stat = 'identity',position = "identity")+
      ggplot2::geom_text(ggplot2::aes(label=Freq),size=3.5)+
      ggplot2::scale_fill_manual(values=pals::coolwarm(2))+
      ggplot2::ggtitle('Receptors')+
      ggplot2::ylab('')+
      ggplot2::xlab('Number of interactions')+
      ggplot2::theme_minimal()
    return(p1+p2)
}


pdf(file = paste0(outdir, 'differential_cci_frequency_table.pdf'), width = 11)
plot_table(all_data)
dev.off()

