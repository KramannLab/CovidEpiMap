# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de


#---- Bar chart of SARS-CoV-2 antibody testing

library(ggplot2)
library(cowplot)
library(dplyr)
library(viridis)
indir = '~/sciebo/CovidEpiMap/antibody_testing/'

data = read.table(file = paste0(indir, 'antibody_test_results.txt'),
                  header = TRUE,
                  sep = '\t')
colnames(data) = c('condition', 'covid_test')
data$condition = factor(data$condition, levels = c('healthy', 'active mild', 'active severe', 
                                                   'recovered mild', 'recovered severe'))

pdf(file = paste0(indir, 'antibody_testing.pdf'), width = 4, height = 4)
ggplot(data, aes(fill = covid_test, x = condition)) + 
  geom_bar(stat = 'count') +
  theme_cowplot() +
  theme(axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size = 8, angle = 90, hjust = 1)) +
  scale_fill_viridis(discrete = TRUE)
dev.off()

