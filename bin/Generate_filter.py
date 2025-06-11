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
	"""
	Parse command line arguments
	"""
	## Standard argparse setup for user input on the command line
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
	## Main function takes the input file and parses it into a format that can be used by OPTIC
	## Essentially, it extracts an filters the required information from the COSMIC cancer mutation census file into a new TSV
	## Checking the input and output directories are formatted correctly so the script can run
	if not arg.output:
		arg.output = os.path.dirname(arg.input) + '/'
	if not os.path.exists(arg.output):
		os.mkdir(arg.output)
	file_path = arg.input
	count = 0
	
	## Want to keep the data as a dictionary so it can be converted to a DataFrame later
	## This is done to avoid having to create a DataFrame and then filter it, which is more memory intensive
	## The dictionary will be converted to a DataFrame at the end of the script
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
	
	## Open file and extract required info. We skip the header line becasue we add a OPTIC specific header later
	with open(file_path, 'rt') as f:
		for line in f:
			count += 1
			if count == 1:
				continue
			else:
				cosmic_information = line.split('\t')
				
				## Skipping lines that do not have the required information becasue it likely indicates improper formatting
				if any(not cosmic_information[i] for i in [0, 1, 2, 3, 6, 7, 15, 16]):
					continue
				# Data to be loaded in final file
				
				gene_id = cosmic_information[0] 
				accession = cosmic_information[1]
				role_in_cancer = cosmic_information[2]
				variant_cds = cosmic_information[6]				
				variant_aa = cosmic_information[7]				
				wt_allele = cosmic_information[11]
				mutation_allele = cosmic_information[12]
				cds_type = cosmic_information[15]
				aa_type = cosmic_information[16]				
				clinvar_sig = cosmic_information[51]				
				cosmic_sig_tier = cosmic_information[57]
				cosmic_tier = str(cosmic_information[3])
				
				## Asigning the "Other" teir to an integer makes filtering with the argparse argument easier
				if cosmic_sig_tier == 'Other\n':
					cosmic_sig_tier = 4
				cosmic_sig_tier = int(cosmic_sig_tier)
				if cosmic_sig_tier > arg.filter_cosmic_sig:
					continue
				
				## Need to ensure the genome build is correct, otherwise the coordinates will be incorrect
				if arg.genome_build == 'GRCh38' or arg.genome_build == 'hg38':
					coordinates = cosmic_information[20]
				else:
					coordinates = cosmic_information[19]
				
				## Splitting the coordinates into chromosome, start and end
				## This is done to ensure the coordinates are in the correct format for OPTIC
				## Skips if the coordinates are not in the correct format"
				result = re.split('[:-]', coordinates)
				if not result:
					continue
				chromosome = result[0]
				start = result[1]
				end = result[2]
				
				# Perform filtering based on the arguments provided
				if arg.filter_t1:
					if cosmic_tier != str(1):
						continue
				if arg.filter_pathogenic:
					if 'pathogenic' not in clinvar_sig.lower():
						continue
				
				## Add to data as dict:
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
	
	## Convert to datarame to save the file
	final_dataframe = pd.DataFrame(data_as_dict)
	final_dataframe.to_csv(f"{arg.output}OPTIC_mutations_file_{arg.genome_build}.tsv", sep="\t", index=False)
	
	
if __name__ == '__main__':
	args = parse_args()

	## Puts the output in the current working directory if no output is specified so files can be found
	if not args.input:
		print('Input folder (--input/-i) is required', file=sys.stderr)
		sys.exit(1)
	if not args.output:
		print('An output location was not provided. Output files will be placed in the input folder')
	
	## Check if the genome build is valid, results will be incorrect or analysis will fail if not
	builds = ['GRCh38', 'hg38', 'GRCh37', 'hg19']
	if args.genome_build not in builds:
		print('The genome build must one of:')
		for item in builds:
			print(f'   {item}')
		sys.exit()
	
	main(args)
	