# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de


#---- TCR richness and evenness across pseudotime 

library(Seurat)
library(viridis)
library(tidyverse)
library(cowplot)
library(vegan)
indir = '~/sciebo/CovidEpiMap/integrated/'
outdir = '~/sciebo/CovidEpiMap/trajectory_analysis/'


# Get pseudotimes per lineage
sc.subset = readRDS(file = paste0(indir, 'integrated.RNA.Tcells.pseudotime.subset.rds'))


# Cut pseudotimes into bins
# Lineage 1
nbins = 8
max.time = max(sc.subset$slingshot_pseudotime_curve1, na.rm = TRUE)
cuts = seq(from = 0, to = max.time, by = max.time / nbins)
sc.subset$slingshot_pseudotime_curve1_cut = cut(sc.subset$slingshot_pseudotime_curve1, breaks = cuts)

# Lineage 2
max.time = max(sc.subset$slingshot_pseudotime_curve2, na.rm = TRUE)
cuts = seq(from = 0, to = max.time, by = max.time / nbins)
sc.subset$slingshot_pseudotime_curve2_cut = cut(sc.subset$slingshot_pseudotime_curve2, breaks = cuts)


# Compute TCR richness and evenness in pseudotime bins (mild and severe)
# Compute per curve, per pseudotime cut, per condition
conditions = c('mild', 'severe')
df = sc.subset@meta.data
df = df[df$condition_collapsed %in% conditions,]
curves = c('slingshot_pseudotime_curve1_cut', 'slingshot_pseudotime_curve2_cut')


data = list()
# Per curve
for (curve in curves){
  # Get cuts for lineage
  cuts = levels(df[,curve])
  tmp.data.list = list()
  
  # Per pseudotime cut
  for (i in 1:length(cuts)){
    cut = cuts[i]
    tmp.list = list()
    
    # Per condition
    for (condition in conditions){
      # Subset to condition and pseudotime cut
      subset = df[df$condition_collapsed %in% condition,]
      # Filter out NA clonotypes
      subset = subset %>% 
        filter(!is.na(TCR_clonotype_id),
               !!sym(curve) == cut)
    
      # Compute relative evenness and richness
      sub.table = table(subset$patient_clonotype)
      evenness = diversity(sub.table, index = 'invsimpson') / length(sub.table)
      richness = length(sub.table) / (sum(sub.table))
      tmp.list[[condition]] = c(curve, cut, evenness, richness, condition)
    }
    
    # Collect score for conditions
    tmp.mat = do.call(rbind, tmp.list)
    rownames(tmp.mat) = NULL
    tmp.data.list[[i]] = tmp.mat
  }
  
  # Collect score for curves
  tmp.data.mat = do.call(rbind, tmp.data.list)
  rownames(tmp.data.mat) = NULL
  data[[curve]] = tmp.data.mat
}


# Transform to data frame
data = as.data.frame(do.call(rbind, data))
colnames(data) = c('curve', 'cut', 'evenness', 'richness', 'condition')
data$evenness = as.numeric(data$evenness)
data$richness = as.numeric(data$richness)

# Get min and max from cut interval
data = data %>% 
      mutate(x_tmp = str_sub(cut, 2, - 2)) %>% 
      separate(x_tmp, c('min', 'max'), sep = ',') %>% 
      mutate_at(c('min', 'max'), as.double)


# Plot TCR evenness
pdf(file = paste0(outdir, 'integrated_Tcells_pseudotime_TCR_evenness.pdf'))
ggplot(data = data, aes(x = max, y = evenness, group = condition)) +
  geom_line(aes(color = condition)) +
  geom_point() +
  scale_colour_viridis(discrete = TRUE, 
                      option = 'viridis',
                      name = 'Condition') + 
  theme_cowplot() + 
  facet_wrap(. ~ curve, 
             labeller = labeller(curve = c(`slingshot_pseudotime_curve1_cut` = 'Lineage 1',
                                `slingshot_pseudotime_curve2_cut` = 'Lineage 2'))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.ticks = element_blank()) +
  labs(x = 'Pseudotime (binned)',
       y = 'TCR evenness')
dev.off()


# Plot TCR richness
pdf(file = paste0(outdir, 'integrated_Tcells_pseudotime_TCR_richness.pdf'))
ggplot(data = data, aes(x = max, y = richness, group = condition)) +
  geom_line(aes(color = condition)) +
  geom_point() +
  scale_colour_viridis(discrete = TRUE, 
                       option = 'viridis',
                       name = 'Condition') + 
  theme_cowplot() + 
  facet_wrap(. ~ curve, 
             labeller = labeller(curve = c(`slingshot_pseudotime_curve1_cut` = 'Lineage 1',
                                           `slingshot_pseudotime_curve2_cut` = 'Lineage 2'))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.ticks = element_blank()) +
  labs(x = 'Pseudotime (binned)',
       y = 'TCR richness')
dev.off()

