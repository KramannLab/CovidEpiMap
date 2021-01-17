# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de


#---- TCR clonal expansion vs pseudotime

library(Seurat)
library(slingshot)
library(ggplot2)
library(viridis)
library(ggbeeswarm)
library(ggthemes)
library(cowplot)
library(gridExtra)
'%ni%' = Negate('%in%')
source('sc_source/sc_source.R')
indir = '~/sciebo/CovidEpiMap/integrated/'
outdir = '~/sciebo/CovidEpiMap/trajectory_analysis/'

# Subset for cells with a clonotype
sc.subset = readRDS(file = paste0(indir, 'integrated.RNA.Tcells.pseudotime.subset.rds'))
sc.subset = subset(sc.subset, patient_clonotype %ni% NA)


sc.subset$clonotype_size_cut = cut(sc.subset$clonotype_size, breaks = c(0, 10, 100, 500, 1000, 1600))
df = sc.subset@meta.data
df = df[df$condition %ni% 'healthy',]
df$integrated_annotations = factor(df$integrated_annotations, 
                                   levels = names(cell.type.colors))


# Lineage 1
# Density
p1 = ggplot() + 
  geom_density(data = df, 
               aes(slingshot_pseudotime_curve1, group = clonotype_size_cut, 
                   fill = clonotype_size_cut), 
               alpha = 0.5, 
               adjust = 2) +
  scale_fill_manual(values = viridis(4), name = 'Clonal expansion') +
  theme_classic() +
  xlab('Pseudotime (Lineage 1)') +
  ylab('Density') +
  theme(axis.ticks = element_blank(),
        axis.text.x = element_blank())

# Cell ordering
p2 = ggplot() +
  geom_point(aes(x = seq_along(df$slingshot_pseudotime_curve1), 
                 y = df$slingshot_pseudotime_curve1, 
                 colour = df$integrated_annotations),
             size = 0.5) +
  scale_colour_manual(values = cell.type.colors) +
  coord_flip() +
  theme_void() +
  NoLegend()


empty_plot = plot(0,type = 'n', axes = FALSE, ann = FALSE)
l = get_legend(p1)
top_row = plot_grid(empty_plot, p2, empty_plot, align = 'h', axis = 'l', ncol = 3, rel_widths = c(1,7.4,2.3))
bottom_row = plot_grid(p1 + NoLegend(), l, align = 'h', axis = 'l', ncol = 2, rel_widths = c(4,1))

pdf(file = paste0(outdir, 'integrated_Tcells_clonal_expansion_lineage1_density.pdf'))
plot_grid(top_row, bottom_row, ncol = 1, rel_heights = c(1/16, 15/16))
dev.off()


# Lineage 2
# Density
p1 = ggplot() + 
  geom_density(data = df, 
               aes(slingshot_pseudotime_curve2, group = clonotype_size_cut, 
                   fill = clonotype_size_cut), 
               alpha = 0.5, 
               adjust = 2) +
  scale_fill_manual(values = viridis(4), name = 'Clonal expansion') +
  theme_classic() +
  xlab('Pseudotime (Lineage 2)') +
  ylab('Density') +
  theme(axis.ticks = element_blank(),
        axis.text.x = element_blank())

# Cell ordering
p2 = ggplot() +
  geom_point(aes(x = seq_along(df$slingshot_pseudotime_curve2), 
                 y = df$slingshot_pseudotime_curve2, 
                 colour = df$integrated_annotations),
             size = 0.5) +
  scale_colour_manual(values = cell.type.colors) +
  coord_flip() +
  theme_void() +
  NoLegend()


l = get_legend(p1)
top_row = plot_grid(empty_plot, p2, empty_plot, align = 'h', axis = 'l', ncol = 3, rel_widths = c(1,7.4,2.3))
bottom_row = plot_grid(p1 + NoLegend(), l, align = 'h', axis = 'l', ncol = 2, rel_widths = c(4,1))

pdf(file = paste0(outdir, 'integrated_Tcells_clonal_expansion_lineage2_density.pdf'))
plot_grid(top_row, bottom_row, ncol = 1, rel_heights = c(1/16, 15/16))
dev.off()



#---- Clonotype size vs pseudotime

df = sc.subset@meta.data
df$integrated_annotations = factor(df$integrated_annotations, 
                                   levels = names(cell.type.colors))

pdf(file = paste0(outdir, 'integrated_Tcells_all_clonotype_size_vs_pseudotime.pdf'), width = 8)
# Lineage 1
ggplot(df, aes(x = slingshot_pseudotime_curve1, y = clonotype_size)) + 
  geom_point(aes(colour = integrated_annotations)) + 
  scale_colour_manual(values = cell.type.colors) +
  stat_smooth(method = 'loess', colour = 'black') +
  facet_wrap(. ~ condition) +
  theme_classic() +
  theme(strip.background = element_rect(colour = 'black', fill = 'lightgrey')) +
  xlab('Pseudotime (Lineage 1)') +
  ylab('Clonotype size') +
  theme(axis.ticks = element_blank(),
        axis.text.x = element_blank())
  
# Lineage 2
ggplot(df, aes(x = slingshot_pseudotime_curve2, y = clonotype_size)) + 
  geom_point(aes(colour = integrated_annotations)) + 
  scale_colour_manual(values = cell.type.colors) +
  stat_smooth(method = 'loess', colour = 'black') +
  facet_wrap(. ~ condition) +
  theme_classic() +
  theme(strip.background = element_rect(colour = 'black', fill = 'lightgrey')) +
  xlab('Pseudotime (Lineage 2)') +
  ylab('Clonotype size') +
  theme(axis.ticks = element_blank(),
        axis.text.x = element_blank())
dev.off()


