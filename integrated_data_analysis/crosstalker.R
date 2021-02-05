# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de


#---- Visualize CellPhoneDB results with CrossTalkeR

source('sc_source/sc_source.R')
source('sc_source/crosstalker_source.R')
source.folder = '~/sciebo/CovidEpiMap/CrossTalkeR/CrossTalkeR/'
devtools::load_all(source.folder)
wkdir = '~/sciebo/CovidEpiMap/CrossTalkeR/'
indir = 'significant_means/filtered_corrected/'
setwd(wkdir)


# Active mild vs active severe
control = 'active_mild'
case = 'active_severe'
data_active = crosstalker_report(Control = control, Case = case, Indir = indir, Outdir = wkdir)


# Recovered mild vs recovered severe
control = 'recovered_mild'
case = 'recovered_severe'
data_recov = crosstalker_report(Control = control, Case = case, Indir = indir, Outdir = wkdir)



#---- Single cell-cell interaction plot (active mild and active severe)

outdir = paste0(wkdir, 'active_mild_vs_active_severe/')

# Active mild
pdf(file = paste0(outdir, 'single_cci_plot_active_mild.pdf'), width = 10, height = 10)
# Modified function to fit coloring scheme
plot_cci(graph = data_active@graphs$CTR,
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

# Active severe
pdf(file = paste0(outdir, 'single_cci_plot_active_severe.pdf'), width = 10, height = 10)
plot_cci(graph = data_active@graphs$EXP,
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



#---- Differential cell-cell interaction plot (active mild vs active severe)

pdf(file = paste0(outdir, 'differential_cci_plot.pdf'), width = 10, height = 10)
plot_cci(graph = data_active@graphs$EXP_x_CTR,
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



#---- Differential cell-cell ligand/receptor frequency table (active mild vs active severe)

all_data = readRDS(file = paste0(outdir, 'LR_data_step2.Rds'))

pdf(file = paste0(outdir, 'differential_cci_frequency_table.pdf'), width = 13)
plot_table(all_data)
dev.off()



#---- Differential gene cell ranking (sending/receiving) (active mild vs active severe)

pdf(file = paste0(outdir, 'differential_gene_cell_ranking.pdf'), width = 14)
plot_rank(all_data)
dev.off()



#----Sankey plots of selected interactions

# Active mild vs active severe
outdir = 'active_mild_vs_active_severe/'
targets = c('MICB', 'IFNG')

for (target in targets){
  pdf(file = paste0(outdir, 'sankey_', target, '.pdf'), width = 10, height = 5)
  print(plot_sankey(data = data_active@tables$EXP_x_CTR, target = target))
  dev.off()
}

# Recovered mild vs recovered severe
outdir = 'recovered_mild_vs_recovered_severe/'
target = 'SELL'

pdf(file = paste0(outdir, 'sankey_', target, '.pdf'), width = 10, height = 5)
plot_sankey(data = data_recov@tables$EXP_x_CTR, target = target)
dev.off()

