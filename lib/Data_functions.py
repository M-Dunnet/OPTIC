"""
Collection of functions to be used for data wrangling in OPTIC.
"""
import sys

import pandas as pd
import datetime
from collections import Counter


def filter_cosmic_mutants(sample_variants_df, cosmic_mutants):
	cosmic_mutants['GENOMIC_WT_ALLELE'] = cosmic_mutants['GENOMIC_WT_ALLELE'].fillna('-')
	cosmic_mutants['GENOMIC_MUT_ALLELE'] = cosmic_mutants['GENOMIC_MUT_ALLELE'].fillna('-')
	
	merged_df = sample_variants_df.merge(cosmic_mutants,
	                                     how='left',
	                                     on='Gene_ID',
	                                     suffixes=('', '_cosmic'))
	merged_df = merged_df.drop(merged_df.index[0])
	merged_df['Position'] = merged_df['Position'].astype(float)
	
	mask = (merged_df['Position'] == merged_df['Position_cosmic']) | (merged_df['Variant'] == merged_df['Variant_cosmic'])
	tmp_df_cosmic = merged_df[mask]
	tmp_df_cosmic.reset_index(drop=True, inplace=True)
	mask = (tmp_df_cosmic['Reference_allele'] == tmp_df_cosmic['GENOMIC_WT_ALLELE']) & (tmp_df_cosmic['Alternate_allele'] == tmp_df_cosmic['GENOMIC_MUT_ALLELE'])
	tmp_df_cosmic = tmp_df_cosmic[mask]
	tmp_df_cosmic = tmp_df_cosmic.drop_duplicates(subset=['Gene_ID', 'Variant'])
	
	return tmp_df_cosmic


def process_optic_dictionary(args, file, sample_variants_df, cosmic_mutants):
	if args.cosmic_mutants:
		sample_variants_df = filter_cosmic_mutants(sample_variants_df, cosmic_mutants)
	
	processed_sample_variants_dict = {}

	for index, row in sample_variants_df.iterrows():
		gene_id = row['Gene_ID']
		variant = row['Variant']
		variant_type = row['Variant_type']
		
		if gene_id not in processed_sample_variants_dict:
			processed_sample_variants_dict[gene_id] = ([variant], [variant_type])
		else:
			processed_sample_variants_dict[gene_id][0].append(variant)
			processed_sample_variants_dict[gene_id][1].append(variant_type)
		
		if gene_id not in processed_sample_variants_dict:
			processed_sample_variants_dict[gene_id] = ([variant], [variant_type])
		else:
			processed_sample_variants_dict[gene_id][0].append(variant)
			processed_sample_variants_dict[gene_id][1].append(variant_type)
	
	gene_dict = {file: processed_sample_variants_dict}
	
	return gene_dict


def count_variant_types(args, variants_dict):
	mutation_types_per_gene = {}
	all_mutations = []
	for filename, gene_data in variants_dict.items():
		for gene, mutation_info in gene_data.items():
			variants = set(mutation_info[0])
			
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
	
	with open(args.output + 'Variant_counts_' + str(datetime.datetime.now().strftime('%Y-%m-%d')) + '.txt', 'w') as variant_log:
		print('Gene' + '\t' + 'Variant' + '\t' + 'Count', file=variant_log)
		for key, value in mutation_counts.items():
			if key == 'Hugo_Symbol':
				continue
			gene, variant = key.split('__')
			count = value
			print(str(gene) + '\t' + str(variant) + '\t' + str(count), file=variant_log)
	
	with open(args.output + 'Variant_types_' + str(datetime.datetime.now().strftime('%Y-%m-%d')) + '.txt', 'w') as variant_types_log:
		print('Gene' + '\t' + 'Variant_Type' + '\t' + 'Count', file=variant_types_log)
		for key, value in mutation_types_per_gene.items():
			if key == 'Hugo_Symbol':
				continue
			gene = key
			for variant_type, count in value.items():
				if variant_type == '':
					variant_type = 'Not Classified'
				print(str(gene) + '\t' + str(variant_type) + '\t' + str(count), file=variant_types_log)


def get_gene_ids(variants_dict):
	target_genes = [gene for value in variants_dict.values() for gene in value.keys()]
	target_genes = set(target_genes)
	if 'Hugo_Symbol' in target_genes:
		target_genes.remove('Hugo_Symbol')
	
	return target_genes


def generate_optic_array(variants_dict, targets):
	gene_list = []
	for gene in targets:
		sample_list = []
		for key, value in variants_dict.items():
			if gene in value:
				sample_list.append(1.0)
			else:
				sample_list.append(0.0)
		gene_list.append(sample_list)
	
	optic_array = pd.DataFrame(gene_list)
	sample_ids = list(variants_dict.keys())
	sample_id_to_column_header = {int(col): sample_ids[col] for col in optic_array.columns}
	optic_array.index = targets
	optic_array['Count'] = optic_array.sum(axis=1)
	optic_array = optic_array.sort_values(by=['Count'], ascending=False)
	
	return optic_array, sample_id_to_column_header


def get_zero_samples(args, mutation_dict):
	zero_number, zero_list = 0, []
	for key in mutation_dict.keys():
		if not mutation_dict[key]:
			zero_number += 1
			zero_list.append(key)
	
	if args.save_zeros:
		zero_list = pd.Series(zero_list)
		zero_list.to_csv(args.output + 'Zero_mutation_samples.tsv', index=False, header=False)
	
	return zero_number


def get_total_frequencies(args, variants_dict, optic_array, maf_list):
	# Initialize Total Gene Mutation Frequencies file
	log_file = open(args.output + 'Gene_mutation_frequencies_' + str(datetime.datetime.now().strftime('%Y-%m-%d')) + '.txt', 'w+')
	zero_number = get_zero_samples(args, variants_dict)
	percentage = round(optic_array['Count'] / len(maf_list) * 100, 2)
	
	print('Total files analysed: ' + str(len(maf_list)), file=log_file)
	print('Files with no observed mutations: ' + str(zero_number), file=log_file)
	print('\n' + 'Gene' + '\t' + 'Mutation Frequency', file=log_file)
	print(percentage.to_string(), file=log_file)
	log_file.close()


def check_optic_array(optic_array):
	if len(optic_array) == 0:
		print(f"Error: OPTIC has not found any valid variants in the input files. Please check that:\n"
		      f"\t (1) Columns numbers are correct in the config file\n "
		      f"\t (2) If a filter file was provided, check that it is suitable for OPTIC and that its for the correct genome build\n")
		sys.exit(1)
		