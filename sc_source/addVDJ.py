#!/usr/bin/env python3
# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de

''' Summarize V, D, J, C genes for each clonotype '''

# Initialise
import sys
clonotypes = dict()
datadir = '/users/monica/Dropbox/CovidEpiMap/data/'
sample = sys.argv[1].upper()
library = sample[sample.find('_')+1:]

wkdir = datadir + sample + '/'
infile = open(wkdir + 'filtered_contig_annotations.csv', 'r')

header = True
for line in infile:
	line = line.strip().split(',')

	# Get positions in header of relevant columns
	if header:
		vpos = line.index('v_gene')
		dpos = line.index('d_gene')
		jpos = line.index('j_gene')
		cpos = line.index('c_gene')
		prodpos = line.index('productive')
		ctypepos = line.index('raw_clonotype_id')
		header = False

	else:
		# Only consider lines with productive clonotypes
		if line[prodpos] == 'True':
			ctype = line[ctypepos]
			vgene = line[vpos]
			dgene = line[dpos]
			jgene = line[jpos]
			cgene = line[cpos]

			# Add new clonotype and genes
			if ctype not in clonotypes:
				clonotypes[ctype] = {library + '_V_GENE': {vgene}, 
				library + '_D_GENE': {dgene}, 
				library + '_J_GENE': {jgene}, 
				library + '_C_GENE': {cgene}}

			# Add new genes to clonotype
			else:
				clonotypes[ctype][library + '_V_GENE'].add(vgene)
				clonotypes[ctype][library + '_D_GENE'].add(dgene)
				clonotypes[ctype][library + '_J_GENE'].add(jgene)
				clonotypes[ctype][library + '_C_GENE'].add(cgene)
infile.close()

# Output
outfile = open(wkdir + 'clonotype.vdj.genes.csv', 'w')
genes = ['V_GENE', 'D_GENE', 'J_GENE', 'C_GENE']
output_genes = [library + '_' + gene for gene in genes]

header_string = 'Clonotype,' + ','.join(output_genes)
outfile.write(header_string + '\n')

for clonotype, genes in clonotypes.items():
	tmp_list = []
	for gene in output_genes:
		tmp_list.append(';'.join(sorted(list(genes[gene]))))

	out_string = clonotype + ',' + ','.join(tmp_list)
	outfile.write(out_string + '\n')
outfile.close()

print('File written:', wkdir + 'clonotype.vdj.genes.csv')


