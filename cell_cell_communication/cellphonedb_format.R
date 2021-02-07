# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de


#---- Format data for input to CellPhoneDB

library(Seurat)
library(tibble)
library(dplyr)
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



#---- Plot CellPhoneDB results


