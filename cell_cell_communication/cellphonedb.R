# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de


#---- Format data for input to CellPhoneDB

library(Seurat)
library(tibble)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(viridis)
indir = '~/sciebo/CovidEpiMap/integrated/'
outdir = '~/sciebo/CovidEpiMap/cellphoneDB/'

sc = readRDS(file = paste0(indir, 'integrated.RNA.Tcells.annotated.rds'))
Idents(sc) = 'condition'
conditions = levels(sc$condition)

for (condition in conditions){
  # Subset data to condition
  subset = subset(sc, condition == condition)
  
  # Format meta data for CellPhoneDB (only keep cell types with count > 20)
  meta = subset@meta.data %>%
    select(integrated_annotations, condition) %>%
    rownames_to_column('sample') %>% 
    group_by(integrated_annotations) %>% 
    add_count %>% 
    filter(n > 20) %>%
    select(-n)
  
  # Fetch normalized count data
  mx = as.matrix(subset@assays$RNA@data)
  
  # Save to file
  write.table(mx, file = paste0(outdir, 'count.', condition, '.txt'), quote = FALSE, sep = '\t')
  write.table(meta, file = paste0(outdir, 'metadata.', condition, '.txt'), quote = FALSE, sep = '\t')
}



# Run cellphonedb on command line
# Example 
# cellphonedb method statistical_analysis metadata.active_mild.txt count.active_mild.txt --counts-data=hgnc_symbol



#---- Plot CellPhoneDB results (From Eleanor Fewings)

condition = 'active_severe'

# Load results
pval = read.table(file = paste(outdir, condition, '/pvalues.txt', sep = ''), 
                  sep = '\t', 
                  stringsAsFactors = FALSE, 
                  header = TRUE, 
                  check.names = FALSE)
means = read.table(file = paste(outdir, condition, '/means.txt', sep = ''), 
                   sep = '\t', 
                   stringsAsFactors = FALSE, 
                   header = TRUE, 
                   check.names = FALSE)
signif = read.table(file = paste(outdir, condition, '/', condition, '_significant_means.txt', sep = ''), 
                    sep = '\t', 
                    stringsAsFactors = FALSE, 
                    header = TRUE, 
                    check.names = FALSE)

# Set rownames
rownames(pval) = pval$id_cp_interaction
rownames(means) = means$id_cp_interaction
rownames(signif) = signif$id_cp_interaction

# Add rank value to means
signif = signif %>% select(id_cp_interaction, rank)
means = merge(signif, means, by = 'id_cp_interaction')

# Create long format
longformat = pval %>% 
  gather(key = 'Pair', 
         value = 'pval', 
         12:ncol(pval)) %>% 
  separate(Pair, c('pair_a', 'pair_b'), 
           sep = '\\|', 
           remove = FALSE)

# Look at interactions between ligand and receptor only
shortened = longformat[(longformat$receptor_a == 'True' & longformat$receptor_b == 'False') | 
                         (longformat$receptor_a == 'False' & longformat$receptor_b == 'True'),]

# Count number of significant interactions
shortened = shortened %>% 
  group_by(Pair) %>% 
  mutate(countint = sum(pval < 0.05))

# Label ligand cell type
shortened$ligand.cell = ifelse(shortened$receptor_a == TRUE, shortened$pair_b, shortened$pair_a)

# Label receptor
shortened$receptor.cell = ifelse(shortened$receptor_a == TRUE, shortened$pair_a, shortened$pair_b)

# Create matrix from data
mx = shortened %>% 
  subset(select = c(ligand.cell, receptor.cell, countint)) %>% 
  unique() %>% 
  mutate(grouped_id = row_number()) %>% 
  spread(receptor.cell, countint) %>% 
  select(-grouped_id)

# Summarise counts
mx = mx %>% 
  group_by(ligand.cell) %>% 
  summarise_all(funs(sum), na.rm = TRUE)

lc = mx$ligand.cell

# Format rownames
mx = mx %>% 
  select(-ligand.cell) %>% 
  as.matrix()

row.names(mx) = lc


# Define colour palette
paletteLength = 100
myColor = colorRampPalette(c(viridis(paletteLength)[1],
                             viridis(paletteLength)[paletteLength]))(paletteLength)

# Plot heatmap
pdf(file = paste0(outdir, condition, '/', condition, '_heatmap_cell_types.pdf'), width = 5, height = 5)
pheatmap(t(mx),
         color = myColor,
         border = NA)
dev.off()

