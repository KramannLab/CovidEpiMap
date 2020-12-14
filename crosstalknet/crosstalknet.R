# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de


#---- Visualize CellPhoneDB with CrossTalkNet (by James Nagai)

# Source scripts
source.folder = '~/sciebo/CovidEpiMap/CrossTalkNet/from_james/CrossTalkNet/'

source.folder = '~/sciebo/Monica_pck/CrossTalkNet/'
devtools::load_all(source.folder) 


# Define path for control and case (CellPhoneDB has been run separately on these)
setwd('~/sciebo/CovidEpiMap/CrossTalkNet/')
indir = 'significant_means/filtered_corrected/'
control = 'active_mild'
case = 'active_severe'

paths = c('CTR' = paste0(indir, control, '_filtered_corrected.csv'),
		'EXP' = paste0(indir, case, '_filtered_corrected.csv'))


# Select gene list
genes = c('MIF', 'KLRK1')


# Generate report
outdir = paste0('/Users/monica/sciebo/CovidEpiMap/CrossTalkNet/', control, '_vs_', case, '/')


data = generate_report(LRpaths = paths,
						genes = genes,
						out_path = outdir,
						threshold = 0,
						out_file = paste0(control, '_vs_', case, '.html'),
						report = TRUE)











