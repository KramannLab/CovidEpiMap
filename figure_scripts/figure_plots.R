# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de


#---- Violin plots of canonical cell type markers

library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(genesorteR)
library(viridis)
'%ni%' = Negate('%in%')
source('sc_source/sc_source.R')


indir = '~/sciebo/CovidEpiMap/integrated/'
sc = readRDS(file = paste0(indir, 'integrated.RNA.Tcells.annotated.rds'))
DefaultAssay(sc) = 'RNA'
Idents(sc) = 'integrated_annotations'


# Genes to plot
genes = read.table(file = paste0(indir, 'marker_expression/selected_cell_type_markers.txt'), header = TRUE)
genes = genes$gene

p = get_violin(object = sc, features.use = genes)
pdf(file = paste0(indir, 'cell_type_markers.pdf'), width = 8, height = 4)
p
dev.off()

get_violin(object = sc, features.use = genes[1:3])



#---- General QC plots

patients = c('29', '31', '32',
             '1', '2', '15',
             '3', '11', '19',
             '5', '6', '7',
             '12', '13', '18')

# Plot cell count per sample in bar chart
data = as.data.frame(table(sc$orig.ident))
colnames(data) = c('sample', 'count')
data$sample = as.character(sub('_GEX_SURF', '', data$sample))
data$sample = factor(data$sample, 
                           levels = patients)

p1 = ggplot(data, aes(x = sample, y = log10(count))) +
  geom_bar(aes(fill = sample), 
           alpha = 1, stat = 'identity') +
  scale_fill_viridis(discrete = TRUE) +
  xlab('') + 
  ylab('Number of valid cells (log10)') +
  theme_cowplot() +
  theme(axis.text.x = element_text(size = 8),
        axis.ticks = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 7),
        legend.position = 'none')

# Plot UMI count per sample in violin plot
data = data.frame(cell = names(Idents(sc)), 
                 UMI = as.numeric(sc$nCount_RNA), 
                 gene = as.numeric(sc$nFeature_RNA),
                 sample = sc$orig.ident)
data$sample = as.character(sub('_GEX_SURF', '', data$sample))
data$sample = factor(data$sample, 
                     levels = patients)

p2 = ggplot(data, aes(x = sample, y = log10(UMI))) +
  geom_boxplot(aes(color = sample),
               outlier.size = 0.5) +
  geom_violin(aes(color = sample, fill = sample), 
              alpha = 0.5) +
  scale_fill_viridis(discrete = TRUE) +
  scale_color_viridis(discrete = TRUE) +
  xlab('') + 
  ylab('Number of UMIs (log10)') +
  theme_cowplot() +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 7),
        axis.ticks = element_blank(),
        legend.position = 'none')

# Plot gene count per sample in violin plot
p3 = ggplot(data, aes(x = sample, y = log10(gene))) +
  geom_boxplot(aes(color = sample),
               outlier.size = 0.5) +
  geom_violin(aes(color = sample, fill = sample), 
              alpha = 0.5) +
  scale_fill_viridis(discrete = TRUE) +
  scale_color_viridis(discrete = TRUE) +
  xlab('') + 
  ylab('Number of genes (log10)') +
  theme_cowplot() +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 7),
        axis.ticks = element_blank(),
        legend.position = 'none')


pdf(file = paste0(indir, 'per_patient_scRNA_QC.pdf'), width = 4, height = 6)
plot_grid(p1, p2, p3, 
          ncol = 1, 
          align = 'v', 
          rel_heights = c(1, 1, 1.3))
dev.off()

#---- Bar charts with average cell type distribution

cell.table = data.frame(cell = colnames(sc), condition = sc$condition,
                        cluster = sc$integrated_annotations, 
                        patient = sc$patient)

df = cell.table %>%
group_by(patient) %>%
count(cluster, name = 'count', .drop = FALSE) %>%
add_count(wt = count, name = 'total') %>%
mutate(fraction = count/total) %>%
as.data.frame

condition = rep(c('healthy', 'active_mild',
        'active_severe', 'recovered_mild',
        'recovered_severe'), 
        each = 3*13)
df = cbind(df, condition)

df = df %>% 
group_by(condition, cluster) %>%
summarize(mean = mean(fraction)) %>%
as.data.frame

df$condition = factor(df$condition, levels = c('healthy', 'active_mild',
        'active_severe', 'recovered_mild',
        'recovered_severe'))
df$condition_collapsed = sub('active_', '', df$condition)
df$condition_collapsed = sub('recovered_', '', df$condition_collapsed)


pdf(file = paste0(indir, 'integrated_Tcells_barchart_celltype_condition_average.pdf'), width = 8)
ggplot(df, aes(fill = cluster, y = mean, x = condition)) + 
    geom_bar(stat = 'identity', position = position_fill(reverse = TRUE)) + 
labs(y = 'Average proportion', x = element_blank(), fill = 'Cell type') +
theme_classic() +
theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
scale_fill_manual(values = cell.type.colors)
dev.off()

pdf(file = paste0(indir, 'integrated_Tcells_barchart_celltype_condition_collapsed_average.pdf'))
ggplot(df, aes(fill = cluster, y = mean, x = condition_collapsed)) + 
    geom_bar(stat = 'identity', position = position_fill(reverse = TRUE)) + 
labs(y = 'Average proportion', x = element_blank(), fill = 'Cell type') +
theme_classic() +
theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
scale_fill_manual(values = cell.type.colors)
dev.off()



#---- Bar chart of CD8+ exhausted T cell distribution

df = cell.table %>% 
  filter(cluster == 'CD8+ exhausted T cells') %>%
  group_by(condition) %>%
  count %>%
  as.data.frame
df = rbind(df, c('active_mild', 0))
df = rbind(df, c('recovered_mild', 0))
df$n = as.integer(df$n)


pdf(file = paste0(indir, 'exhausted_tcells_count.pdf'), width = 4, height = 3.5)
ggplot(df, aes(x = condition, y = n, fill = condition)) +
  geom_bar(stat = 'identity') +
  theme_cowplot() +
  scale_fill_viridis(discrete = TRUE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.ticks = element_blank(),
        axis.title.x = element_blank()) +
  labs(y = 'Cell type count', 
       fill = 'Condition') +
  scale_y_continuous(breaks = c(0,25, 50))
dev.off()



#---- Heatmap of ADT library expression

DefaultAssay(sc) = 'ADT'
markers = rownames(sc)[24:35]

pdf(file = paste0(indir, 'integrated_Tcells_ADT_selected_expression.pdf'), height = 5.5, width = 8)
DoHeatmap(subset(sc, downsample = 500), 
          features = markers, 
          assay = 'ADT', 
          raster = FALSE, 
          group.colors = cell.type.colors,
          group.bar.height = 0.03, 
          label = FALSE) +
  scale_fill_gradient2(low = viridis(2)[1], 
                       mid = '#F0F0F0', 
                       high = viridis(2)[2], 
                       limits = c(-1,1), 
                       oob = scales::squish)
dev.off()


