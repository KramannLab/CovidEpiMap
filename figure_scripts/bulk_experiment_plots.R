# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de


#---- Tile plot of bulk seq experiment, mode of dextramer binding
library(ggplot2)
library(ggsci)
library(dplyr)
library(viridis)
library(reshape2)
library(magrittr)
library(writexl)
library(tidyverse)
library(cowplot)
library(scales)

indir = '~/sciebo/CovidEpiMap/bulk_experiment/'
setwd(indir)


# Get dextramer binding specificity and HLA compatibility
Match_enrichment = read.table('comparison_HLA_alleles_patient_vs_dextramer.txt', header = TRUE, sep = '\t')

Match_enrichment_compare = mutate(Match_enrichment,
                                  mode.of.dextramer.binding = ifelse(as.character(match) == TRUE & log2_Enrichment >= 1, 'Specific enrichment',
                                                              ifelse(as.character(match) == TRUE & log2_Enrichment < 1, 'No antigen specificity',
                                                              ifelse(as.character(match) == FALSE & log2_Enrichment >= 1, 'Unspecific binding',
                                                              ifelse(as.character(match) == FALSE & log2_Enrichment < 1, 'Not HLA compatible', NA)))))

# Format columns
Match_enrichment_compare$peptide.name = paste(Match_enrichment_compare$HLA.allele,
                                              Match_enrichment_compare$epitope.peptide.sequence)
Match_enrichment_compare$sample.name = sub('SCV2-', 'P', Match_enrichment_compare$Sample)

# Plot
pdf(file = 'dextramer_enrichment_tile.pdf', width = 13, height = 8)
ggplot(Match_enrichment_compare) +
  geom_tile(mapping = aes(x = sample.name, y = sub('\\d+ ', '', peptide.name), fill = mode.of.dextramer.binding), color = 'black') + 
  facet_grid(rows = vars(HLA.supertype), cols = vars(infection.status), scales = 'free', space = 'free') +
  scale_fill_viridis(discrete = TRUE, option = 'viridis') +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        legend.title = element_text(size = 8),
        strip.text.x = element_text(size = 8),
        axis.text.x = element_text(color = 'black', size = 8, angle = 90, vjust = 0.5),
        axis.text.y = element_text(color = 'black', size = 8),
        legend.text = element_text(size = 8),
        axis.ticks = element_blank()) +
  labs(fill = 'Dextramer binding')
dev.off()



#---- Bar chart of bulk seq experiment, dextramer enrichment

# Compute mean dextramer enrichment
specific_enrichment = Match_enrichment_compare %>% 
  filter(mode.of.dextramer.binding == 'Specific enrichment') %>%
  group_by(peptide.name, infection.status) %>%
  summarize(mean = mean(log2_Enrichment),
            sd = sd(log2_Enrichment))
specific_enrichment$infection.status = as.factor(specific_enrichment$infection.status)

# Get HLA supertype
match.supertype = unique(Match_enrichment_compare %>% select(peptide.name, HLA.supertype))
hla.supertype = match.supertype$HLA.supertype
names(hla.supertype) = match.supertype$peptide.name
specific_enrichment$hla.supertype = hla.supertype[specific_enrichment$peptide.name]
specific_enrichment$hla.supertype = factor(specific_enrichment$hla.supertype, levels = c('A02', 'A03', 'A01'))
specific_enrichment$infection.status= factor(specific_enrichment$infection.status, 
                                             levels = c('healthy', 'mild active', 'severe active',
                                                        'mild recovered', 'severe recovered'))


# Plot
pdf(file = 'dextramer_enrichment_bar.pdf',width = 8, height = 4)
ggplot(specific_enrichment, aes(x = reorder(sub('\\d+ ', '', peptide.name), -mean, 'mean'), 
                                y = mean, 
                                fill = infection.status,
                                ymin = mean - sd,
                                ymax = mean + sd)) + 
  geom_bar(stat = 'identity', position = position_dodge2(width = 0.9, preserve = 'single')) +
  geom_errorbar(stat = 'identity',
                position = position_dodge2(width = 0.9, preserve = 'single')) +
  scale_fill_viridis(discrete = TRUE, option = 'viridis') +
  theme_classic() +
  facet_grid(~ hla.supertype, scales = 'free', space = 'free', switch = 'x') +
  theme(axis.text.x = element_text(color = 'black', size = 8, angle = 90, vjust = 0.5),
        axis.ticks = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  xlab('SARS-COV2 epitope') +
  ylab(bquote(Log[2]~'(enrichment)')) +
  labs(fill = 'Condition')
dev.off()

