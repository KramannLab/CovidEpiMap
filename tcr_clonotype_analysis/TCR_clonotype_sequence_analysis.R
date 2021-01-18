# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de


#---- TCR sequence analysis of A0101-2 binding cells

library(Seurat)
library(ggplot2)
library(tidyr)
library(dplyr)
library(tidytext)
library(viridis)
library(msa)
library(seqinr)
indir = '~/sciebo/CovidEpiMap/integrated/'
outdir = '~/sciebo/CovidEpiMap/tcr_sequence_analysis/'
  
dex.subset = readRDS(file = paste0(indir, 'A0101-2.unique.binders.rds'))
df = dex.subset@meta.data
df = df[,c('condition_collapsed', 'patient_clonotype', 'TCR_cdr3s_aa', 'TCR_V_GENE', 'TCR_J_GENE')]

# Count CDR3 sequences of epitope-binding cells
data = df %>% 
  group_by(condition_collapsed) %>% 
  dplyr::count(TCR_cdr3s_aa, name = 'TCR_cdr3s_aa_count') %>% 
  arrange(desc(TCR_cdr3s_aa_count)) %>% 
  as.data.frame()

write.table(data,
            file = paste0(outdir, 'A0101-2.unique.binders.cdr3s.counts.txt'),
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



#---- Multiple alignments of top15 TRA and TRB sequences across conditions

# Save top 15 TRA and TRB sequences in fasta format
tcr.seq.data = tcr.data %>% 
  group_by(TCR_chain, condition_collapsed) %>% 
  top_n(n = 15, wt = percent) %>% 
  as.data.frame()

chains = c('TRA', 'TRB')
conditions = c('healthy', 'mild', 'severe')

for (chain in chains){
  # Subset for CR chain
  tmp.data = tcr.seq.data %>%
    filter(TCR_chain == chain)
  
  # Subset for condition 
  for (condition in conditions){
    tmp = tmp.data %>%
      filter(condition_collapsed == condition)
    
    # Save sequences as fasta files
    tcr.seq = tmp$TCR_cdr3s_aa
    tcr.names = paste0(condition, 1:15)
    datdir = paste0(outdir, chain,  '_fasta/')
    
    for (i in 1:15){
      write.fasta(sequences = tcr.seq[i],
                  names = tcr.names[i],
                  file.out = paste0(datdir, tcr.names[i], '.fasta'))
    }
  }
}


# Multiple alignment with Clustal Omega
setwd(paste0(outdir, 'clustal_omega/'))

for (chain in chains){
  for (condition in conditions){
    # Read sequences
    datdir = paste0(outdir, chain,  '_fasta/')
    fasta = paste0(datdir, condition, 1:15, '.fasta')
    sequences = readAAStringSet(fasta)
    
    # Alignment
    alignment = msaClustalOmega(inputSeqs = sequences, 
                                type = 'protein', 
                                order = 'input')
    
    # Output alignments and consensus to pdf
    msaPrettyPrint(alignment,
                   output = 'pdf',
                   file = paste0(chain, '_', condition, '_multiple_alignment.pdf'),
                   showNames = 'left',
                   showLogo = 'top',
                   shadingMode = 'similar',
                   askForOverwrite = FALSE,
                   furtherCode = c('\\defconsensus{.}{lower}{upper}',
                                   '\\showruler{1}{top}'))
  }
}



#---- Investigate changes in TCR usage over pseudotime

# Add pseudotimes to A0101-2 binding cells
sc.subset = readRDS(file = paste0(indir, 'integrated.RNA.Tcells.pseudotime.subset.rds'))

df = sc.subset@meta.data
dex.subset.cells = colnames(dex.subset)
dex.subset$slingshot_pseudotime = df[dex.subset.cells, 'slingshot_pseudotime']
dex.subset$slingshot_pseudotime_curve1 = df[dex.subset.cells, 'slingshot_pseudotime_curve1']
dex.subset$slingshot_pseudotime_curve2 = df[dex.subset.cells, 'slingshot_pseudotime_curve2']


# Cut pseudotimes into bins
nbins = 8
max.time = max(dex.subset$slingshot_pseudotime, na.rm = TRUE)
cuts = seq(from = 0, to = max.time, by = max.time / nbins)
dex.subset$slingshot_pseudotime_cut = cut(dex.subset$slingshot_pseudotime, breaks = cuts)

cut.names = levels(cut(dex.subset$slingshot_pseudotime, breaks = cuts))
names(cut.names) = paste0('pseudotime_cut_', 1:nbins)
dex.subset$slingshot_pseudotime_cut_name = names(cut.names[dex.subset$slingshot_pseudotime_cut])
dex.subset = subset(dex.subset, slingshot_pseudotime_cut %ni% NA)


# Get TRA and TRB sequences
df = dex.subset@meta.data
df = df %>% 
  separate_rows(TCR_cdr3s_aa, sep = ';') %>%
  as.data.frame
df$TCR_chain = 'TRA'
df[grep('TRB:', df$TCR_cdr3s_aa),'TCR_chain'] = 'TRB'
df$TCR_cdr3s_aa = gsub('TR(A|B):', '', df$TCR_cdr3s_aa)


tcr.data = df %>% 
  group_by(slingshot_pseudotime_cut_name, TCR_chain) %>% 
  dplyr::count(TCR_cdr3s_aa) %>% 
  arrange(desc(n)) %>%
  top_n(n = 15, wt = n) %>%
  as.data.frame


# Save top 15 TRA and TRB sequences in fasta format
for (chain in chains){
  # Subset for CR chain
  tmp.data = tcr.data %>%
    filter(TCR_chain == chain)
  
  # Subset for pseudotime bin 
  for (cut in names(cut.names)){
    tmp = tmp.data %>%
      filter(slingshot_pseudotime_cut_name == cut)
    
    # Save sequences as fasta files
    tcr.seq = tmp$TCR_cdr3s_aa
    tcr.names = paste0(cut, '_', 1:15)
    datdir = paste0(outdir, chain,  '_fasta/')
    
    for (i in 1:15){
      write.fasta(sequences = tcr.seq[i],
                  names = tcr.names[i],
                  file.out = paste0(datdir, tcr.names[i], '.fasta'))
    }
  }
}


# Multiple alignment with Clustal Omega
setwd(paste0(outdir, 'clustal_omega/'))

for (chain in chains){
  for (cut in names(cut.names)){
    # Read sequences
    datdir = paste0(outdir, chain,  '_fasta/')
    fasta = paste0(datdir, cut, '_', 1:15, '.fasta')
    sequences = readAAStringSet(fasta)
    
    # Alignment
    alignment = msaClustalOmega(inputSeqs = sequences, 
                                type = 'protein', 
                                order = 'input')
    
    # Output alignments and consensus to pdf
    msaPrettyPrint(alignment,
                   output = 'pdf',
                   file = paste0(chain, '_', cut, '_multiple_alignment.pdf'),
                   showNames = 'left',
                   showLogo = 'top',
                   shadingMode = 'similar',
                   askForOverwrite = FALSE,
                   furtherCode = c('\\defconsensus{.}{lower}{upper}',
                                   '\\showruler{1}{top}'))
  }
}



#---- Visualise pseudotime bins on UMAP
pdf(paste0(outdir, 'A0101-2_unique_binders_pseudotime_cuts.pdf'), width = 8.7)
DimPlot(dex.subset, group.by = 'slingshot_pseudotime_cut_name', cols = viridis(nbins)) + 
  labs(colour = 'Pseudotime bin') +
  theme(axis.ticks = element_blank())
dev.off()


