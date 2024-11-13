#!/usr/bin/env python3
"""
The script can return functions from 'gene_panel_functions.py'.
"""

import os
import sys
import pandas as pd
import argparse
import datetime
from collections import Counter

PATH = '/'.join(os.path.realpath(__file__).split('/')[:-1]) + '/lib/'
sys.path.append(os.path.abspath(PATH))

import gene_panel_functions as optic


def parse_args():
	parser = argparse.ArgumentParser(description='Assesses the ability of a panel to cover all samples.',
	                                 add_help=True,
	                                 prefix_chars='-')
	# Required options
	parser.add_argument('--input', '-i', type=str, action='store',
	                    help='The directory location containing input maf files')
	parser.add_argument('--output', '-o', type=str, action='store',
	                    help='The directory location where files will end up')
	
	# Analysis Options
	parser.add_argument('--targets', '-t', type=str, action='store',
	                    help='Path string to the location of the target file. Only loci in the target file will be '
	                         'analysed')
	parser.add_argument('--greedy_coverage', action='store_true', default=False,
	                    help='Boolean option to run the greedy set coverage algorythm. This will select the subsets of a given gene '
	                         'to cover as many samples as possible, prioritizing subsets that cover the most uncovered elements at each step.'
	                         'Only recommended for targeted panels')
	parser.add_argument('--hierarchial_clustering', type=int, action='store', default=None,
	                    help='Sets oncoscan to cluster samples by SNP and INDEL profiles. The integer provided sets the number of clusters to analyse')
	
	# Filter options
	parser.add_argument('--cosmic_mutants', type=str, action='store',
	                    help='Path string to the location of the cosmic mutant census file. Only mutations in the cosmic mutation census'
	                         ' will be analysed')
	parser.add_argument('--numGenes', '-n', type=int, action='store', default=25,
	                    help='Maximum number of genes to examine for the mutation matrix, set-cover, and joint occurrence matrix. Defaults to 25, maximum of 50')
	parser.add_argument('--numClusterGenes', '-c', type=int, action='store', default=250,
	                    help='Maximum number of genes to examine during hierarchial clustering. Defaults to 250. Note this number cannot be lower than numGenes.')
	parser.add_argument('--save_zeros', '-z', default=False, action='store_true',
	                    help='Boolean option to save a list of all samples that do not contain mutations in the gene set')
	parser.add_argument('--use_cds', default=False, action='store_true',
	                    help='If set, Oncoscan will use CDS changes rather than protein changes')
	
	# Plotting Options
	parser.add_argument('--fixed_matrix_width', '-f', default=False, action='store_true',
	                    help='Boolean option to plot the mutation matrix where each sample has a fixed width. '
	                         'If not set, column widths will be scaled relative to the sample number. Default is False')
	parser.add_argument('--plot_total_variants', default=False, action='store_true',
	                    help='Boolean option to plot the total number of variants per sample (includes passenger mutations). Default is False')
	return parser.parse_args()


def targeted_oncoscan(arg, variants_dict, targets, maf_list, log_file):
	if arg.cosmic_mutants:
		mutation_types_per_gene = {}
		all_mutations = []
		for filename, gene_data in variants_dict.items():
			for gene, mutation_info in gene_data.items():
				variants = mutation_info[0]
				for cds_variant in variants:
					mutations = str(gene) + '__' + str(cds_variant)
					all_mutations.append(mutations)
				
				mutation_types = mutation_info[1]
				counter = Counter(mutation_types)
				if gene not in mutation_types_per_gene:
					mutation_types_per_gene[gene] = counter
				else:
					mutation_types_per_gene[gene] += counter
		mutation_counts = Counter(all_mutations)
		
		with open(arg.output + 'Variant_counts_' + str(datetime.datetime.now().strftime('%Y-%m-%d')) + '.txt', 'w') as cds_variant_log:
			print('Gene' + '\t' + 'Variant' + '\t' + 'Count', file=cds_variant_log)
			for key, value in mutation_counts.items():
				gene, variant = key.split('__')
				count = value
				print(str(gene) + '\t' + str(variant) + '\t' + str(count), file=cds_variant_log)
		
		with open(arg.output + 'Variant_types_' + str(datetime.datetime.now().strftime('%Y-%m-%d')) + '.txt', 'w') as variant_types_log:
			print('Gene' + '\t' + 'Variant_Type' + '\t' + 'Count', file=variant_types_log)
			for key, value in mutation_types_per_gene.items():
				gene = key
				for variant_type, count in value.items():
					if variant_type == '':
						variant_type = 'Not Classified'
					print(str(gene) + '\t' + str(variant_type) + '\t' + str(count), file=variant_types_log)
	
	gene_list = []
	for gene in targets:
		sample_list = []
		for key, value in variants_dict.items():
			if gene in value:
				sample_list.append(1.0)
			else:
				sample_list.append(0.0)
		gene_list.append(sample_list)
	df = pd.DataFrame(gene_list)
	sample_ids = list(variants_dict.keys())
	sample_id_to_column_header = {int(col): sample_ids[col] for col in df.columns}
	df.index = targets
	
	df['Count'] = df.sum(axis=1)
	df = df.sort_values(by=['Count'], ascending=False)
	percentage = round(df['Count'] / len(maf_list) * 100, 2)
	df = df.iloc[:, :-1]
	
	print('\n' + 'Gene' + '\t' + 'Mutation Frequency', file=log_file)
	print(percentage.to_string(), file=log_file)
	log_file.close()
	
	return df, sample_id_to_column_header


def all_genes_oncoscan(arg, variants_dict, maf_list, log_file):
	if arg.cosmic_mutants:
		mutation_types_per_gene = {}
		all_mutations = []
		
		for filename, gene_data in variants_dict.items():
			for gene, mutation_info in gene_data.items():
				variants = mutation_info[0]
				for variant in variants:
					mutations = str(gene) + '__' + str(variant)
					all_mutations.append(mutations)
				
				mutation_types = mutation_info[1]
				counter = Counter(mutation_types)
				if gene not in mutation_types_per_gene:
					mutation_types_per_gene[gene] = counter
				else:
					mutation_types_per_gene[gene] += counter
		mutation_counts = Counter(all_mutations)
		
		with open(arg.output + 'Variant_counts_' + str(datetime.datetime.now().strftime('%Y-%m-%d')) + '.txt', 'w') as variant_log:
			print('Gene' + '\t' + 'Variant' + '\t' + 'Count', file=variant_log)
			for key, value in mutation_counts.items():
				gene, variant = key.split('__')
				count = value
				print(str(gene) + '\t' + str(variant) + '\t' + str(count), file=variant_log)
		
		with open(arg.output + 'Variant_types_' + str(datetime.datetime.now().strftime('%Y-%m-%d')) + '.txt', 'w') as variant_types_log:
			print('Gene' + '\t' + 'Variant_Type' + '\t' + 'Count', file=variant_types_log)
			for key, value in mutation_types_per_gene.items():
				gene = key
				for variant_type, count in value.items():
					if variant_type == '':
						variant_type = 'Not Classified'
					print(str(gene) + '\t' + str(variant_type) + '\t' + str(count), file=variant_types_log)
	
	all_genes_list = [gene for value in variants_dict.values() for gene in value.keys()]
	all_genes_list = set(all_genes_list)
	if 'Hugo_Symbol' in all_genes_list:
		all_genes_list.remove('Hugo_Symbol')
	
	gene_list = []
	for gene in all_genes_list:
		sample_list = []
		for key, value in variants_dict.items():
			if gene in value:
				sample_list.append(1.0)
			else:
				sample_list.append(0.0)
		gene_list.append(sample_list)
	df = pd.DataFrame(gene_list)
	sample_ids = list(variants_dict.keys())
	sample_id_to_column_header = {int(col): sample_ids[col] for col in df.columns}
	df.index = all_genes_list
	
	df['Count'] = df.sum(axis=1)
	df = df.sort_values(by=['Count'], ascending=False)
	percentage = round(df['Count'] / len(maf_list) * 100, 2)
	df = df.iloc[:, :-1]
	
	print(percentage.to_string(), file=log_file)
	log_file.close()
	
	if arg.numGenes > 50:
		print('The Maximum number of genes has been set too high! --numGenes/-n has now been set to 50.')
	
	return df, sample_id_to_column_header


def main(arg):
	maf_list = []
	os.chdir(arg.input)
	date = str(datetime.datetime.now().strftime('%Y-%m-%d'))
	if arg.targets:
		arg.output = f'{arg.output}Targeted_oncoscan_{date}/'
		if not os.path.exists(arg.output):
			os.mkdir(arg.output)
	elif arg.hierarchial_clustering:
		arg.output = f'{arg.output}Hierarchial_clustering_{arg.hierarchial_clustering}_clusters_{date}/'
		if not os.path.exists(arg.output):
			os.mkdir(arg.output)
	else:
		arg.output = f'{arg.output}All_genes_oncoscan_{date}/'
		if not os.path.exists(arg.output):
			os.mkdir(arg.output)
	
	# -------------------------------------------------------------------------------------
	# Change as necessary based on file extension.
	#  -------------------------------------------------------------------------------------
	for maf_file in os.listdir(arg.input):
		if maf_file.endswith('.tsv'):
			maf_list.append(maf_file)
	
	maf_list = list(set(maf_list))
	maf_list.sort()
	
	if len(maf_list) == 0:
		print('There are no VAF files in input folder...')
		sys.exit(1)
	
	log_file = open(arg.output + 'Gene_mutation_frequencies_' + str(datetime.datetime.now().strftime('%Y-%m-%d')) + '.txt', 'w+')
	
	if args.targets:
		with open(arg.targets, 'rt') as f:
			targets = f.read().splitlines()
	
	# Get list of mutated genes
	mutation_dict, number = {}, 0
	for file in maf_list:
		number += 1
		raw_sample_variants_dict = {}
		processed_sample_variants_dict = {}
		
		with open(file, 'rt') as f:
			for i, line in enumerate(f):
				if not line.startswith('#'):
					columns = line.strip().split('\t')
					# -------------------------------------------------------------------------------------
					# Change as necessary based VCF file format. Removes variants based on specific column content
					#  -------------------------------------------------------------------------------------
					if columns[0] == 'Unknown':
						continue
					if columns[10] == 'synonymous_variant':
						continue
					if arg.use_cds:
						raw_sample_variants_dict[i] = [columns[0], columns[5], columns[11], columns[13], columns[34]]
					else:
						raw_sample_variants_dict[i] = [columns[0], columns[5], columns[11], columns[13], columns[39]]
					# columns[0] = Gene,
					# columns[5] = Position,
					# columns[11] = Reference allele,
					# columns[13] = Alternate allele,
					# columns[34] = Variant CDS (only used when --use_cds is applied)
					# columns[39] = variant amino acid change
		
		if arg.cosmic_mutants:
			cosmic_mutants = pd.read_csv(arg.cosmic_mutants, sep='\t')
			cosmic_mutants = cosmic_mutants.astype({'Position': float, 'Variant_CDS': str, 'Variant_AA': str})
			if arg.use_cds:
				cosmic_mutants.rename(columns={'Variant_CDS': 'Variant'}, inplace=True)
			else:
				cosmic_mutants.rename(columns={'Variant_AA': 'Variant'}, inplace=True)
				
			# -------------------------------------------------------------------------------------
			# Change as necessary based on INDEL reference/allele format of the VCF files
			# (i.e. deletions might be recorded as an empty string (e.g. '') or as a hyphen ('-').
			#  -------------------------------------------------------------------------------------
			cosmic_mutants['GENOMIC_WT_ALLELE'] = cosmic_mutants['GENOMIC_WT_ALLELE'].fillna('-')
			cosmic_mutants['GENOMIC_MUT_ALLELE'] = cosmic_mutants['GENOMIC_MUT_ALLELE'].fillna('-')

			sample_variants_df = pd.DataFrame.from_dict(raw_sample_variants_dict, orient='index')
			sample_variants_df.columns = ['Gene_ID', 'Position', 'Reference_allele', 'Alternate_allele', 'Variant']
			
			merged_df = sample_variants_df.merge(cosmic_mutants,
			                                     how='left',
			                                     on='Gene_ID',
			                                     suffixes=('_sample', '_cosmic'))
			merged_df = merged_df.drop(merged_df.index[0])
			merged_df['Position_sample'] = merged_df['Position_sample'].astype(float)
			
			mask = (merged_df['Position_sample'] == merged_df['Position_cosmic']) | (merged_df['Variant_sample'] == merged_df['Variant_cosmic'])
			tmp_df_cosmic = merged_df[mask]
			tmp_df_cosmic.reset_index(drop=True, inplace=True)
			mask = (tmp_df_cosmic['Reference_allele'] == tmp_df_cosmic['GENOMIC_WT_ALLELE']) & (tmp_df_cosmic['Alternate_allele'] == tmp_df_cosmic['GENOMIC_MUT_ALLELE'])
			tmp_df_cosmic = tmp_df_cosmic[mask]
			tmp_df_cosmic = tmp_df_cosmic.drop_duplicates(subset=['Gene_ID', 'Variant_cosmic'])

			for index, row in tmp_df_cosmic.iterrows():
				gene_id = row['Gene_ID']
				variant = row['Variant_cosmic']
				if arg.use_cds:
					variant_type = row['CDS_Type']
				else:
					variant_type = row['AA_Type']
				if not isinstance(variant_type, str):
					variant_type = ''
				role_in_cancer = row['Role_in_cancer']
				
				if gene_id not in processed_sample_variants_dict:
					processed_sample_variants_dict[gene_id] = ([variant], [variant_type], role_in_cancer)
				else:
					processed_sample_variants_dict[gene_id][0].append(variant)
					processed_sample_variants_dict[gene_id][1].append(variant_type)
			
			if arg.targets:
				processed_sample_variants_dict = {key: value for key, value in processed_sample_variants_dict.items() if key in targets}
				gene_dict = {file: processed_sample_variants_dict}  # switch file/number
			else:
				gene_dict = {file: processed_sample_variants_dict}
			
			mutation_dict.update(gene_dict)
		
		else:
			sample_variants_df = pd.DataFrame.from_dict(raw_sample_variants_dict, orient='index')
			sample_variants_df.columns = ['Gene_ID', 'Position', 'Reference_allele', 'Alternate_allele', 'Variant']
			for index, row in sample_variants_df.iterrows():
				gene_id = row['Gene_ID']
				variant = row['Variant']
				
				if gene_id not in processed_sample_variants_dict:
					processed_sample_variants_dict[gene_id] = ([variant])
				else:
					processed_sample_variants_dict[gene_id].append(variant)
			if arg.targets:
				processed_sample_variants_dict = {key: value for key, value in processed_sample_variants_dict.items() if key in targets}
				gene_dict = {file: processed_sample_variants_dict}  # switch file/number
			else:
				gene_dict = {file: processed_sample_variants_dict}
			mutation_dict.update(gene_dict)
		
		print(f"Finished file {number} of {len(maf_list)} ({round(number / len(maf_list) * 100, 2)}%)", ' ' * 20, end='\r')
	
	zero_number, zero_list = 0, []
	for key in mutation_dict.keys():
		if not mutation_dict[key]:
			zero_number += 1
			zero_list.append(key)
	print('Total files analysed: ' + str(len(maf_list)), file=log_file)
	print('Files with no observed mutations: ' + str(zero_number), file=log_file)
	if arg.save_zeros:
		zero_list = pd.Series(zero_list)
		zero_list.to_csv(arg.output + 'Zero_mutation_samples.tsv', index=False, header=False)
	
	if arg.targets:
		df, sample_mapping = targeted_oncoscan(arg, mutation_dict, targets, maf_list, log_file)
	else:
		df, sample_mapping = all_genes_oncoscan(arg, mutation_dict, maf_list, log_file)
	
	gsc_df = optic.greedy_set_cover(df, arg.numGenes)
	
	if args.hierarchial_clustering:
		optic.hierarchial_clustering(arg, df, sample_mapping, args.hierarchial_clustering)
	
	elif args.greedy_coverage:
		optic.clusterplot(arg, gsc_df)
		optic.gene_set_coverage(arg, gsc_df)
		optic.joint_occurrence_matrix(arg, gsc_df)
		optic.plot_column_sums(arg, gsc_df)
	
	else:
		optic.clusterplot(arg, df)
		optic.gene_set_coverage(arg, df)
		optic.joint_occurrence_matrix(arg, df)
		optic.plot_column_sums(arg, df)
	
	if args.plot_total_variants:
		optic.total_sample_mutations(arg, maf_list)


if __name__ == '__main__':
	args = parse_args()
	if not args.input:
		print('Error: an input folder (--input/-i) is required', file=sys.stderr)
		sys.exit(1)
	if not args.output:
		print('Error: an output folder (--output/-o) is required', file=sys.stderr)
		sys.exit(1)
	if not args.output.endswith('/'):
		args.output = args.output + '/'
	if args.targets and args.hierarchial_clustering:
		print('Error: A targeted list of genes (--targets/-t) and hierarchial clustering (--hierarchial_clustering) cannot be provided simultaneously.')
		sys.exit(1)
	if args.greedy_coverage and args.hierarchial_clustering:
		print('Error: The greedy set coverage (--greedy_coverage) algorithm and hierarchial clustering (--hierarchial_clustering) cannot be provided simultaneously.')
		sys.exit(1)
	if args.greedy_coverage and not args.targets:
		print('Warning: The greedy set coverage algorithm (--greedy_coverage) is designed for checking targeted panels. It is unlikely to perform well on all genes')
	main(args)
