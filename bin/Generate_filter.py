#!/usr/bin/env python3
"""
Formats the COSMIC cancer mutation census file into a readable version for OncoScan
"""

import os
import sys
import pandas as pd
import argparse
import re


def parse_args():
	parser = argparse.ArgumentParser(description='Formats the COSMIC cancer mutation census file into a readable version for OncoScan.',
	                                 add_help=True,
	                                 prefix_chars='-')
	parser.add_argument('--input', '-i', type=str, action='store',
	                    help='File path of the cancer mutation census tsv file')
	parser.add_argument('--output', '-o', type=str, action='store',
	                    help='Location to save parsed DataFrame')
	parser.add_argument('--genome_build', '-g', type=str, action='store',
	                    help='Defines which genome version to use, must be either `GRCh38` (or `hg38`) or `GRCh37` (or `hg19`)')
	parser.add_argument('--filter_t1', '-t', action='store_true', default=False,
	                    help='Boolean. If set, all genes not within the tier 1 COSMIC gene census will be excluded. Defaults to False')
	parser.add_argument('--filter_pathogenic', '-p', action='store_true', default=False,
	                    help='Boolean. If set, only pathogenic ClinVar mutations will be kept. Defaults to False')
	parser.add_argument('--filter_cosmic_sig', '-c', action='store', type=int, default=4,
	                    help='Set as INT between 1-4. If set, only mutations at that tier level and higher will be kept (e.g. if set as 2, tier 1 and 2 mutations will be kept.'
	                         'The `Other` tier is renamed to 4. default is 4 (No filtering)')
	
	return parser.parse_args()


def main(arg):
	if not arg.output:
		arg.output = os.path.dirname(arg.input) + '/'
	if not os.path.exists(arg.output):
		os.mkdir(arg.output)
	# set file location of the COSMIC cancer mutation census file
	file_path = arg.input
	count = 0
	
	data_as_dict = {'Gene_ID': [],
	                'Accession': [],
	                'Chromosome': [],
	                'Position': [],
	                'CDS_Type': [],
	                'AA_Type': [],
	                'Variant_CDS': [],
	                'Variant_AA': [],
	                'GENOMIC_WT_ALLELE': [],
	                'GENOMIC_MUT_ALLELE': [],
	                'Role_in_cancer': [],
	                'ClinVar_Sig': [],
	                'Cosmic_Sig': []}
	
	# Open file and extract required info
	with open(file_path, 'rt') as f:
		for line in f:
			count += 1
			if count == 1:
				continue
			else:
				cosmic_information = line.split('\t')
				
				# Data to be loaded in final file
				gene_id = cosmic_information[0]
				if not gene_id:
					continue
				accession = cosmic_information[1]
				if not accession:
					continue
				role_in_cancer = cosmic_information[2]
				if not role_in_cancer:
					continue
				variant_cds = cosmic_information[6]
				if not variant_cds:
					continue
				variant_aa = cosmic_information[7]
				if not variant_aa:
					continue
				wt_allele = cosmic_information[11]
				mutation_allele = cosmic_information[12]
				cds_type = cosmic_information[15]
				if not cds_type:
					continue
				aa_type = cosmic_information[16]
				if not aa_type:
					continue
				clinvar_sig = cosmic_information[51]

				cosmic_sig_tier = cosmic_information[57]
				if cosmic_sig_tier == 'Other\n':
					cosmic_sig_tier = 4
				cosmic_sig_tier = int(cosmic_sig_tier)
				if cosmic_sig_tier > arg.filter_cosmic_sig:
					continue
				
				if arg.genome_build == 'GRCh38' or arg.genome_build == 'hg38':
					coordinates = cosmic_information[20]
				else:
					coordinates = cosmic_information[19]
				result = re.split('[:-]', coordinates)
				if not result:
					continue
				chromosome = result[0]
				start = result[1]
				end = result[2]

				# Other data for filtering
				cosmic_tier = str(cosmic_information[3])
				if not cosmic_tier:
					continue
				
				# Filtering
				if arg.filter_t1:
					if cosmic_tier != str(1):
						continue
				
				if arg.filter_pathogenic:
					if 'pathogenic' not in clinvar_sig.lower():
						continue
				
				# Add to data as dict:
				data_as_dict['Gene_ID'].append(gene_id)
				data_as_dict['Accession'].append(accession)
				data_as_dict['Chromosome'].append(chromosome)
				data_as_dict['Position'].append(start)
				data_as_dict['CDS_Type'].append(cds_type)
				data_as_dict['AA_Type'].append(aa_type)
				data_as_dict['Variant_CDS'].append(variant_cds)
				data_as_dict['Variant_AA'].append(variant_aa)
				data_as_dict['GENOMIC_WT_ALLELE'].append(wt_allele)
				data_as_dict['GENOMIC_MUT_ALLELE'].append(mutation_allele)
				data_as_dict['Role_in_cancer'].append(role_in_cancer)
				data_as_dict['ClinVar_Sig'].append(clinvar_sig)
				data_as_dict['Cosmic_Sig'].append(cosmic_sig_tier)
	
	final_dataframe = pd.DataFrame(data_as_dict)
	final_dataframe.to_csv(f"{arg.output}OPTIC_mutations_file_{arg.genome_build}.tsv", sep="\t", index=False)
	
	
if __name__ == '__main__':
	args = parse_args()
	if not args.input:
		print('Input folder (--input/-i) is required', file=sys.stderr)
		sys.exit(1)
	if not args.output:
		print('An output location was not provided. Output files will be placed in the input folder')
	
	builds = ['GRCh38', 'hg38', 'GRCh37', 'hg19']
	if args.genome_build not in builds:
		print('The genome build must one of:')
		for item in builds:
			print(f'   {item}')
		sys.exit()
	
	main(args)
	