# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de


#---- Visualize CellPhoneDB results with CrossTalkNet (by James Nagai)

# Source scripts
source.folder = '~/sciebo/CovidEpiMap/CrossTalkNet/from_james/CrossTalkNet/'
devtools::load_all(source.folder)

wkdir = '~/sciebo/CovidEpiMap/CrossTalkNet/'
indir = 'significant_means/filtered_corrected/'
setwd(wkdir)


crosstalknet_report = function(Control, Case, Genes, Indir, Outdir){
	# Define path for control and case (CellPhoneDB has been run separately on these)
	# Files have already been filtered and corrected with organize.py from CrossTalkNet package
	paths = c('CTR' = paste0(Indir, Control, '_filtered_corrected.csv'),
			'EXP' = paste0(Indir, Case, '_filtered_corrected.csv'))

	# Generate report
	outdir = paste0(Outdir, '/', Control, '_vs_', Case, '/')

	generate_report(LRpaths = paths,
					genes = Genes,
					out_path = outdir,
					threshold = 0,
					out_file = paste0(Control, '_vs_', Case, '.html'),
					report = TRUE)
}


# Healthy vs active mild
control = 'healthy'
case = 'active_mild'
genes = c('MIF', 'KLRK1')
crosstalknet_report(Control = control, Case = case, Genes = genes, Indir = indir, Outdir = wkdir)


# Healthy vs active severe
control = 'healthy'
case = 'active_severe'
genes = c('MIF', 'KLRK1')
crosstalknet_report(Control = control, Case = case, Genes = genes, Indir = indir, Outdir = wkdir)


# Healthy vs recovered mild
control = 'healthy'
case = 'recovered_mild'
genes = c('MIF', 'KLRK1')
crosstalknet_report(Control = control, Case = case, Genes = genes, Indir = indir, Outdir = wkdir)


# Healthy vs recovered severe
control = 'healthy'
case = 'recovered_severe'
genes = c('MIF', 'KLRK1')
crosstalknet_report(Control = control, Case = case, Genes = genes, Indir = indir, Outdir = wkdir)


# Active mild vs active severe
control = 'active_mild'
case = 'active_severe'
genes = c('MIF', 'KLRK1')
crosstalknet_report(Control = control, Case = case, Genes = genes, Indir = indir, Outdir = wkdir)


