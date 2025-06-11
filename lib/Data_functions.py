"""
Collection of functions to be used for data wrangling in OPTIC.
"""
import sys
import pandas as pd
import datetime
from collections import Counter


def filter_cosmic_mutants(sample_variants_df, cosmic_mutants):
	"""
	Filter the sample variants DataFrame to only include variants that match those in the COSMIC mutants DataFrame.
	"""
	## Data is in Pandas DataFrame format, so we can use the merge function to filter
	## NA in the COSMIC DataFrame is being filled with '-' so that it matches the sample_variants_df
	cosmic_mutants['GENOMIC_WT_ALLELE'] = cosmic_mutants['GENOMIC_WT_ALLELE'].fillna('-')
	cosmic_mutants['GENOMIC_MUT_ALLELE'] = cosmic_mutants['GENOMIC_MUT_ALLELE'].fillna('-')
	
	## Merge on Gene ID; this leads to multiple mismatched rows, but these will be filtered out later.
	## Also remove the header row from the merged DataFrame; after merge we end up with a duplicate header row
	## Poisition is converted to float to ensure that it can be compared correctly
	merged_df = sample_variants_df.merge(cosmic_mutants,
	                                     how='left',
	                                     on='Gene_ID',
	                                     suffixes=('', '_cosmic'))
	merged_df = merged_df.drop(merged_df.index[0])
	merged_df['Position'] = merged_df['Position'].astype(float)
	
	## Filter the merged DataFrame to only include rows where the Position and Variant match those in the COSMIC DataFrame
	## This removes mismatched rows where the Position or Variant do not match
	mask = (merged_df['Position'] == merged_df['Position_cosmic']) | (merged_df['Variant'] == merged_df['Variant_cosmic'])
	tmp_df_cosmic = merged_df[mask]
	tmp_df_cosmic.reset_index(drop=True, inplace=True)

	## Filter the DataFrame to only include rows where the Reference and Alternate alleles match those in the COSMIC DataFrame
	## This further removes mismatched rows where the Reference or Alternate alleles do not match; but position and variant do
	mask = (tmp_df_cosmic['Reference_allele'] == tmp_df_cosmic['GENOMIC_WT_ALLELE']) & (tmp_df_cosmic['Alternate_allele'] == tmp_df_cosmic['GENOMIC_MUT_ALLELE'])
	tmp_df_cosmic = tmp_df_cosmic[mask]
	tmp_df_cosmic = tmp_df_cosmic.drop_duplicates(subset=['Gene_ID', 'Variant'])
	
	return tmp_df_cosmic


def process_optic_dictionary(args, file, sample_variants_df, cosmic_mutants):
	"""
	Process the sample variants DataFrame to create a dictionary of gene IDs and their associated variants and types.
	"""
	## Filter first, if a filter file was provided
	if args.cosmic_mutants:
		sample_variants_df = filter_cosmic_mutants(sample_variants_df, cosmic_mutants)
	
	## Formats the DataFrame to have Gene_ID, Variant, and Variant_type columns
	## We collect this information for the Variant counts and Variant types files 
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
	
	## Gene dict is returned, for later use in a nested dictionary format
	gene_dict = {file: processed_sample_variants_dict}
	
	return gene_dict


def count_variant_types(args, variants_dict):
	"""
	Processes the variants dictionary (Nested dictionary made from process_optic_dictionary) to count the number of variants and their types per gene, across all samples.
	"""
	## Iterates through each sample, and collects variant counts and variant types
	## Collect these separately, so that we can write them to different files later
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
	
	## Write the counts to files
	with open(args.output + 'Variant_counts_' + str(datetime.datetime.now().strftime('%Y-%m-%d')) + '.txt', 'w') as variant_log:
		print('Gene' + '\t' + 'Variant' + '\t' + 'Count', file=variant_log)
		for key, value in mutation_counts.items():
			if key == 'Hugo_Symbol': ## Exlclude Hugo_Symbol if it is present; this is a header row that is not a gene
				continue
			gene, variant = key.split('__')
			count = value
			print(str(gene) + '\t' + str(variant) + '\t' + str(count), file=variant_log)
	
	with open(args.output + 'Variant_types_' + str(datetime.datetime.now().strftime('%Y-%m-%d')) + '.txt', 'w') as variant_types_log:
		print('Gene' + '\t' + 'Variant_Type' + '\t' + 'Count', file=variant_types_log)
		for key, value in mutation_types_per_gene.items():
			if key == 'Hugo_Symbol': ## Exlclude Hugo_Symbol if it is present; this is a header row that is not a gene
				continue
			gene = key
			for variant_type, count in value.items():
				if variant_type == '':
					variant_type = 'Not Classified'
				print(str(gene) + '\t' + str(variant_type) + '\t' + str(count), file=variant_types_log)


def get_gene_ids(variants_dict):
	"""
	Collects all unique gene IDs from the variants dictionary.
	"""
	target_genes = [gene for value in variants_dict.values() for gene in value.keys()]
	target_genes = set(target_genes)
	if 'Hugo_Symbol' in target_genes: ## Exlclude Hugo_Symbol if it is present; this is a header row that is not a gene
		target_genes.remove('Hugo_Symbol')
	
	return target_genes


def generate_optic_array(variants_dict, targets):
	"""
	Generates a binary mutation matrix (OPTIC array) from the variants dictionary.
	Each row corresponds to a gene, and each column corresponds to a sample.
	1 indicates the presence of a mutation in that gene for that sample, 0 indicates absence.
	"""
	## For each gene, measure if it is mutated in each sample
	gene_list = []
	for gene in targets:
		sample_list = []
		for key, value in variants_dict.items():
			if gene in value:
				sample_list.append(1.0)
			else:
				sample_list.append(0.0)
		gene_list.append(sample_list)
	
	## Convert the list of lists into a DataFrame
	## Dataframe is used because it makes plotting easier later on
	optic_array = pd.DataFrame(gene_list)
	sample_ids = list(variants_dict.keys())
	sample_id_to_column_header = {int(col): sample_ids[col] for col in optic_array.columns}
	optic_array.index = targets
	optic_array['Count'] = optic_array.sum(axis=1) ## Sorted by the number of mutations per gene
	optic_array = optic_array.sort_values(by=['Count'], ascending=False)
	
	return optic_array, sample_id_to_column_header


def get_zero_samples(args, mutation_dict):
	"""
	Counts the number of samples with no mutations and saves the list of these samples if requested."""
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
	"""
	Calculates and logs the total mutation frequencies across all samples and genes.
	"""
	log_file = open(args.output + 'Gene_mutation_frequencies_' + str(datetime.datetime.now().strftime('%Y-%m-%d')) + '.txt', 'w+')
	zero_number = get_zero_samples(args, variants_dict)
	percentage = round(optic_array['Count'] / len(maf_list) * 100, 2)
	
	print('Total files analysed: ' + str(len(maf_list)), file=log_file)
	print('Files with no observed mutations: ' + str(zero_number), file=log_file)
	print('\n' + 'Gene' + '\t' + 'Mutation Frequency', file=log_file)
	print(percentage.to_string(), file=log_file)
	log_file.close()


def check_optic_array(optic_array):
	"""
	Checks if the OPTIC array is empty, which indicates that no valid variants were found.
	"""
	if len(optic_array) == 0:
		print(f"Error: OPTIC has not found any valid variants in the input files. Please check that:\n"
		      f"\t (1) Columns numbers are correct in the config file\n "
		      f"\t (2) If a filter file was provided, check that it is suitable for OPTIC and that its for the correct genome build\n")
		sys.exit(1)
		