import os
import sys
import pandas as pd
import argparse


def parse_args():
	parser = argparse.ArgumentParser(description='Limits the mutation filter file to only contain variants from specific regions at specified genes.',
	                                 add_help=True,
	                                 prefix_chars='-')
	parser.add_argument('--input', '-i', type=str, action='store',
	                    help='The location of the Parsed COSMIC mutants tsv file')
	parser.add_argument('--hotspots', type=str, action='store',
	                    help='The location of the hotspots file. Must be in TSV format')
	parser.add_argument('--output', '-o', type=str, action='store',
	                    help='Location to save hotspot filtered dataframe')
	return parser.parse_args()


def main(args):
	cosmic_mutants_path = args.input
	hotspots_path = args.hotspots
	
	parsed_cosmic_mutants = pd.read_csv(cosmic_mutants_path, sep='\t', dtype={'Gene_ID': str,
	                                                                          'Accession': str,
	                                                                          'Chromosome': str,
	                                                                          'position': float,
	                                                                          'CDS_Type': str,
	                                                                          'AA_Type': str,
	                                                                          'Variant_CDS': str,
	                                                                          'Variant_AA': str,
	                                                                          'GENOMIC_WT_ALLELE': str,
	                                                                          'GENOMIC_MUT_ALLELE': str,
	                                                                          'Role_in_cancer': str,
	                                                                          'ClinVar_Sig': str,
	                                                                          'Cosmic_Sig': int})
	hotspots_dataframe = pd.read_csv(hotspots_path, sep='\t')
	
	merged = pd.merge(parsed_cosmic_mutants, hotspots_dataframe, on='Gene_ID', how='left', suffixes=('', '_df2'))

	filtered_hotspots = merged[(merged['Position'] >= merged['Start']) & (merged['Position'] <= merged['End'])]

	hotspot_gene_set = set(hotspots_dataframe['Gene_ID'])
	mask = ~parsed_cosmic_mutants['Gene_ID'].isin(hotspot_gene_set)
	result = parsed_cosmic_mutants[mask]
	
	result = pd.concat([result, filtered_hotspots], axis=0)
	result = result.reset_index(drop=True)
	print(result)
	result = result.iloc[:, :13]
	print(result)
	result.fillna('', inplace=True)
	print(result)
	result.to_csv(args.output + 'OncoScan_mutations_file_Hotspots.tsv', sep="\t", index=False)
	

if __name__ == '__main__':
	args = parse_args()
	if not args.input:
		print('Input file (--input/-i) is required', file=sys.stderr)
		sys.exit(1)
	if not args.hotspots:
		print('hotspots file (--hotspots/-s) is required', file=sys.stderr)
		sys.exit(1)
	if not args.output:
		print('Output folder (--output/-o) is required', file=sys.stderr)
		sys.exit(1)
	if not args.output.endswith('/'):
		args.output = args.output + '/'
	
	main(args)
