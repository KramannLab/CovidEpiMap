# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de


#---- TCR sequence analysis of A0101-2 binding cells

library(Seurat)
library(ggplot2)
library(tidyr)
library(dplyr)
library(tidytext)
library(viridis)
indir = '~/sciebo/CovidEpiMap/integrated/'
outdir = '~/sciebo/CovidEpiMap/tcr_sequence_analysis/'
  
dex.subset = readRDS(file = paste0(indir, 'A0101-2.unique.binders.rds'))
df = dex.subset@meta.data

# Count CDR3 sequences of epitope-binding cells (per patient)
data = df %>% 
  group_by(patient) %>% 
  dplyr::count(TCR_cdr3s_aa, name = 'TCR_cdr3s_aa_count') %>% 
  arrange(desc(TCR_cdr3s_aa_count)) %>% 
  as.data.frame()

write.table(data,
            file = paste0(outdir, 'A0101-2.unique.binders.cdr3s.counts.patient.txt'),
            sep = '\t',
            row.names = FALSE,
            quote = FALSE)

# Count CDR3 sequences of epitope-binding cells (per condition collapsed)
data = df %>% 
  group_by(condition_collapsed) %>% 
  dplyr::count(TCR_cdr3s_aa, name = 'TCR_cdr3s_aa_count') %>% 
  arrange(desc(TCR_cdr3s_aa_count)) %>% 
  as.data.frame()

write.table(data,
            file = paste0(outdir, 'A0101-2.unique.binders.cdr3s.counts.condition.collapsed.txt'),
            sep = '\t',
            row.names = FALSE,
            quote = FALSE)


# Plot TRA and TRB frequencies per condition
df = df %>% 
  separate_rows(TCR_cdr3s_aa, sep = ';') %>%
  as.data.frame
df$TCR_chain = 'TRA'
df[grep('TRB:', df$TCR_cdr3s_aa),'TCR_chain'] = 'TRB'
df$TCR_cdr3s_aa = gsub('TR(A|B):', '', df$TCR_cdr3s_aa)

# Summarise data
tcr.data = df %>% 
  group_by(condition_collapsed, TCR_chain) %>% 
  dplyr::count(TCR_cdr3s_aa) %>% 
  arrange(desc(n)) %>% 
  group_by(TCR_chain) %>%
  mutate(percent = (n / sum(n))*100)


# Plot
pdf(file = paste0(outdir, 'A0101-2_unique_binders_TRB_seq_frequency.pdf'), height = 5)
tcr.data %>% 
  filter(TCR_chain == 'TRB') %>% 
  group_by(condition_collapsed) %>% 
  top_n(n = 15, wt = percent) %>%
  ggplot(aes(x = reorder_within(TCR_cdr3s_aa, rev(percent), condition_collapsed), 
             y = percent, 
             fill = condition_collapsed)) +
  geom_bar(stat = 'identity', show.legend = FALSE) +
  scale_x_reordered() +
  cowplot::theme_cowplot() +
  theme(axis.ticks = element_blank(),
        axis.text.x = element_text(size = 8, 
                                   angle = 90, 
                                   hjust = 1, 
                                   vjust = 0.5),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_blank()) +
  facet_grid(~ condition_collapsed, 
             scales = 'free_x', 
             space = 'free') +
  scale_fill_viridis(discrete = TRUE, option = 'viridis') +
  labs(y = 'TRB (%)')
dev.off()

