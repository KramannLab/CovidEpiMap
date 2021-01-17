# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de


#---- Sequence similarity between dextramer epitopes and other COV vira

library(msa)
library(seqinr)
library(viridis)
library(reshape2)
outdir = '~/sciebo/CovidEpiMap/epitope_sequence_analysis/'

# Read dextramer names and sequences
dextramers = read.table(file = '~/sciebo/CovidEpiMap/bulk_experiment/comparison_HLA_alleles_patient_vs_dextramer.txt',
                        sep = '\t', header = TRUE)
dextramers = unique(dextramers[,c('HLA.allele', 'epitope.peptide.sequence')])

# Format name
dextramers$HLA.allele = gsub(' $', '', dextramers$HLA.allele)
dextramers$HLA.allele = gsub(' ', '_', dextramers$HLA.allele)
dextramers$HLA.allele = gsub('\\*', '', dextramers$HLA.allele)
dextramers$HLA.allele = paste0('dex', dextramers$HLA.allele)

# Write out as fasta format
dex.names = dextramers$HLA.allele
dex.seq = dextramers$epitope.peptide.sequence
datdir = paste0(outdir, 'dextramer_fasta/')

for (i in 1:length(dex.names)){
  write.fasta(sequences = dex.seq[i], 
              names = dex.names[i], 
              file.out = paste0(datdir, dex.names[i], '.fasta'))
}


# Multiple alignment with ClustalOmega

viral.fasta = list.files(paste0(outdir, 'viral_proteome_fasta'), full.names = TRUE)
dex.fasta = list.files(datdir, full.names = TRUE)
# Do not include positive CMV controls
dex.fasta = dex.fasta[1:34]
# Do not include dextramers not used in scRNA experiment
dex.fasta = dex.fasta[c(1:6, 9:11, 17:19, 21, 27:28, 31:34)]
dex.names = dex.names[c(1:6, 9:11, 17:19, 21, 27:28, 31:34)]

# Add protein source of dextramer epitope sequences
dex.source = c('Spike', 'Spike', 'Membrane', 'Nucleoprotein',
               'Nucleoprotein', 'Nucleoprotein', 'Spike', 'Spike',
               'Spike', 'Spike', 'Spike', 'Spike',
               'Membrane', 'Spike', 'Spike', 'Negative', 
               'Negative', 'Negative', 'Negative')


setwd(paste0(outdir, 'clustal_omega/'))
for (i in 1:length(dex.fasta)){
  # Get dextramer and viral sequences
  sequences = readAAStringSet(c(dex.fasta[i], viral.fasta))
  
  # Subset sequences to relevant protein source
  if (dex.source[i] ==  'Negative'){
    # Check against spike protein for negative control
    sequences = sequences[c(1, grep('Spike', names(sequences)))]
  } else {
    sequences = sequences[c(1, grep(dex.source[i], names(sequences)))]
  }
  names(sequences) = gsub("\\s*\\[[^\\)]+\\]","", names(sequences))
  
  # Multiple alignment with Clustal Omega
  alignment = msaClustalOmega(inputSeqs = sequences, type = 'protein', order = 'input')
  
  
  # Output alignments and consensus to pdf
  msaPrettyPrint(alignment,
                 output = 'pdf',
                 file = paste0(dex.names[i], '_multiple_alignment.pdf'),
                 showNames = 'left',
                 showLogo = 'none',
                 shadingMode = 'similar',
                 askForOverwrite = FALSE,
                 furtherCode = c('\\defconsensus{.}{lower}{upper}',
                                 '\\showruler{1}{top}'))
  
  
  # Compute similarity (1 - distance) based on alignment results (Fitch matrix)
  alignment.res = msaConvert(alignment, type = 'seqinr::alignment')
  sim.mat = 1 - dist.alignment(alignment.res, 'similarity')
  sim.mat = as.matrix(sim.mat)
  coef_matrix = melt(get_upper_triangle(sim.mat, diag = FALSE), na.rm = TRUE)
  
  # Plot similarity
  pdf(file = paste0(outdir, 'dextramer_similarity/', dex.names[i], '_similarity.pdf'))
  print(upper_trig_tile_plot(coef_matrix, text_size = 3.5) +
    labs(fill = 'Protein similarity'))
  dev.off()
}

