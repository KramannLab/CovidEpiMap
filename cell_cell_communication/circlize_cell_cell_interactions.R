# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de


#---- Plot selected ligand-receptor interactions from CellPhoneDB/CrossTalkeR (single active conditions)

library(circlize)
library(dplyr)
library(viridis)
library(ggplot2)
library(ComplexHeatmap)
source('sc_source/sc_source.R')
indir = '~/sciebo/CovidEpiMap/CrossTalkeR/active_mild_vs_active_severe/'
outdir = '~/sciebo/CovidEpiMap/cell_cell_interaction/'
condition = 'active_severe'


# Define selected interactions
interactions = c('CD94:NKG2A', 'CD94:NKG2C', 
                 'NKG2D II receptor', 'CD94:NKG2E',
                 'KIR3DL2','HLA-C', 'HLA-C')

names(interactions) = c('HLA-E', 'HLA-E',
                        'MICB', 'HLA-E',
                        'HLA-F','KIR2DL1', 'KIR2DL3')

# Get abbreviated cell type names for plotting (remove TEX)
cell.types = names(cell.type.colors[-11])
names(cell.types) = c('TN', 'TCM', 'Treg',
                      'TEMRA', 'NK TEMRA', 'TEM1',
                      'TEM2', 'Tcyc', 'NK TEFF',
                      'NKTatyp', 'MAIT', 'GDT')


# Read data
data = read.csv(file = paste0(indir, 'single_', condition, '.csv'), header = TRUE)
data$X = NULL
data$Ligand.Cluster = names(cell.types[match(data$Ligand.Cluster, cell.types)])
data$Receptor.Cluster = names(cell.types[match(data$Receptor.Cluster, cell.types)])
data$from_type = paste(data$Ligand.Cluster, 'L', sep = ':')
data$to_type = paste(data$Receptor.Cluster, 'R', sep = ':')
colnames(data) = c('ligand.cluster', 'receptor.cluster', 'from', 'to', 'value', 'from_type', 'to_type')
data$value = as.numeric(sub(',', '.', data$value))


# Subset selected ligand-receptor pairs
sub = paste(data$to, data$from, sep = ':') %in% paste(interactions, names(interactions), sep = ':')
subset = data[sub,]


# Get number of outer sections (cell types)
sectors = names(cell.types)
nsectors = length(sectors)


# Define  width for cell types in outer track based on total meanlr values
meanl = subset %>% dplyr::group_by(ligand.cluster) %>% summarize(meanl = sum(value))
meanr = subset %>% dplyr::group_by(receptor.cluster) %>% summarize(meanr = sum(value))
ml = meanl$meanl
names(ml) = meanl$ligand.cluster
mr = meanr$meanr
names(mr) = meanr$receptor.cluster
mlr = c(ml, mr)
width = tapply(mlr, names(mlr), sum)


# Get number of inner sections (ligand/receptor) and order 
#table(c(subset$from_type, subset$to_type))
lr.sectors = c('TN:L', 'TN:R', 
               'TCM:L', 'TCM:R',
               'Treg:L', 'Treg:R',
               'TEMRA:L', 'TEMRA:R', 
               'NK TEMRA:L', 'NK TEMRA:R', 
               'TEM1:L', 'TEM1:R', 
               'TEM2:L', 'TEM2:R', 
               'Tcyc:L', 'Tcyc:R',
               'NK TEFF:L', 'NK TEFF:R',
               'NKTatyp:L', 'NKTatyp:R', 
               'MAIT:L', 'MAIT:R',
               'GDT:L', 'GDT:R')

n.lr.sectors = length(lr.sectors)

# Define widths for inner track
names(ml) = paste0(names(ml), ':L')
names(mr) = paste0(names(mr), ':R')
lr.width = c(ml, mr)


# Define colors for outer track
# Exclude TEX
cell.type.colors = cell.type.colors[-11]
names(cell.type.colors) = names(cell.types)

# Define colors for inner track
lr.sector.colors = ifelse(grepl(':R', lr.sectors, fixed = TRUE), viridis(3)[1], viridis(3)[2])
names(lr.sector.colors) = lr.sectors  

# Define colors for links
yellow = rgb(255, 255, 51, alpha = 200, max = 255) # HLA-E > CD94:NKG2C
darkred = rgb(139, 0, 0, alpha = 200, max = 255) # MICB > NKG2D II receptor
pink = rgb(255, 182, 193, alpha = 200, max = 255) # HLA-E > CD94:NKG2E
lavender = rgb(202, 185, 241, alpha = 200, max = 255) # HLA-E > CD94:NKG2A
darkblue = rgb(0, 0, 88, alpha = 200, max = 255) # HLA-F > KIR3DL2
purple = rgb(128, 0, 128, alpha = 200, max = 255) # KIR2DL1 > HLA-C
indigo = rgb(75, 0, 130, alpha = 200, max = 255) # KIR2DL3 > HLA-C

interaction.colors = c(yellow, darkred, pink, lavender, darkblue, purple, indigo)
names(interaction.colors) = c('HLA-E:CD94:NKG2C',
                              'MICB:NKG2D II receptor',
                              'HLA-E:CD94:NKG2E',
                              'HLA-E:CD94:NKG2A',
                              'HLA-F:KIR3DL2',
                              'KIR2DL1:HLA-C',
                              'KIR2DL3:HLA-C')

# Function to plot outer and inner track
tracks = function(){
  # Outer track
  circos.clear()
  gap = 5
  minigap = 0.7
  circos.par(gap.degree = gap)
  circos.initialize(sectors = sectors, xlim = cbind(rep(0, nsectors), width))
  circos.track(sectors = sectors, 
               ylim = c(0,1),
               track.height = mm_h(8), 
               bg.col = cell.type.colors,
               bg.border = NA)
  
  # Inner track
  circos.clear()
  par(new = TRUE)
  gaps = rep(c(minigap, gap), n.lr.sectors/2)
  circos.par('canvas.xlim' = c(-1.15, 1.15), 'canvas.ylim' = c(-1.15, 1.15), gap.degree = gaps,
             cell.padding = c(0.02, 0, 0.02, 0))
  circos.initialize(sectors = lr.sectors, xlim = cbind(rep(0, n.lr.sectors), lr.width))
  circos.track(sectors = lr.sectors, 
               ylim = c(0,1),
               track.height = mm_h(8),
               bg.col = lr.sector.colors,
               bg.border = NA)
}


# Define plot legends
lr.legend = Legend(at = c('Ligand', 'Receptor'),
                   type = 'grid', 
                   legend_gp = gpar(fill = c(viridis(3)[2], viridis(3)[1])),
                   title_position = 'topleft', 
                   title = 'Inner Track')

link.legend = Legend(at = c('HLA-E - CD94:NKG2C', 'MICB - NKG2D II receptor',
                            'HLA-E - CD94:NKG2E', 'HLA-E - CD94:NKG2A',
                            'HLA-F - KIR3DL2', 'HLA-C - KIR2DL1',
                            'HLA-C - KIR2DL3'),
                     type = 'grid',
                     legend_gp = gpar(fill = c(yellow, darkred,
                                               pink, lavender,
                                               darkblue, purple,
                                               indigo)),
                     title_position = 'topleft',
                     title = 'Ligand - Receptor interaction')

packed.legends = packLegend(lr.legend, link.legend)


# Initiate list for start and end position of receptor links
link.list = list()
for (lr in names(lr.width)){
  link.list[[lr]] = c('start' = unname(lr.width[lr]), 'end' = unname(lr.width[lr]))
}

# Keep track of wheter it's the first link to the receptor
first.lr = rep(TRUE, length(lr.width))
names(first.lr) = names(lr.width)



# Plot
# Add outer and inner tracks
pdf(file = paste0(outdir, 'single_', condition, '_interactions.pdf'), width = 10, height = 10)
tracks()
# Loop over ligands to add links
ligands = names(lr.width)[grepl(':L', names(lr.width))]
for (ligand in ligands){
  # Subset interactions for ligand
  end1 = 0
  dat = subset %>% filter(from_type == ligand) %>% arrange(desc(to_type))
  for (i in 1:nrow(dat)){
    # Get color for interaction
    int = paste(dat[i,'from'], dat[i,'to'], sep = ':')
    int.col = interaction.colors[int]
    
    # Start and end position for ligand cluster
    start1 = end1
    end1 = end1 + dat[i, 'value']
    from = dat[i,'from_type']
    to = dat[i,'to_type']
    
    # Start and end position for receptor cluster 
    if(first.lr[to]){
      start2 = link.list[[to]]['start']
      end2 = link.list[[to]]['end'] - dat[i, 'value']
      circos.link(from, c(start1, end1), to, c(start2, end2), col = int.col)
      link.list[[to]] = c(start2 - dat[i, 'value'], end2)
      first.lr[to] = FALSE
    } else{
      start2 = link.list[[to]]['start']
      end2 = link.list[[to]]['end'] - dat[i, 'value']
      circos.link(from, c(start1, end1), to, c(start2, end2), col = int.col)
      link.list[[to]] = c(start2 - dat[i, 'value'], end2)
    }
  }
}

# Print legends
draw(packed.legends, x = unit(6, 'mm'), y = unit(4, 'mm'), just = c('left', 'bottom'))
dev.off()

