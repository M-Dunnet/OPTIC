"""
Collection of analysis functions to be used in OPTIC.
"""
import os
import pandas as pd
import numpy as np
import datetime
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, linkage


def gene_set_coverage(args, data_matrix, output_modifier=''):
	"""
    Creates a line plot of the number of genes in a panel vs total set coverage (%).
    Uses a greedy set coverage algorithm (from function: greedy_set_cover) .
    """
	if len(data_matrix.index) > args.numGenes:
		data_matrix = data_matrix.head(args.numGenes)
	
	gene_number, zero_number, gene_ids = [], [], []
	gene_number.append(int(0))
	zero_number.append(int(0))
	gene_ids.append(' ')
	for index in data_matrix.index:
		gene_ids.append(index)
	
	iteration = range(1, (data_matrix.shape[0]) + 1)
	for i in iteration:
		dat = data_matrix.head(i)
		num = dat.sum().eq(0).sum()
		data_point = (1 - (num / data_matrix.shape[1])) * 100
		zero_number.append(data_point)
		gene_number.append(i)
	
	set_coverage_data = list(zip(gene_number, zero_number))
	set_coverage_dataframe = pd.DataFrame(set_coverage_data, columns=['Gene_ID', 'Cumulative_percent_coverage'])
	set_coverage_dataframe.to_csv(args.output + output_modifier + 'Set_coverage_datafile_' + str(datetime.datetime.now().strftime('%Y-%m-%d')) + '.tsv', sep='\t', index=False)
	
	max_font_size = 16
	font_size = min(max_font_size, 500 / len(data_matrix))
	
	plt.figure(figsize=(8, 6))
	sns.set_style('whitegrid', {'grid.color': 'lightgrey'})
	sns.lineplot(
		x=gene_number,
		y=zero_number,
		color='black')
	plt.xlabel('Cumulative Gene Number', fontsize=16)
	plt.xlim(0, len(gene_number))
	plt.xticks(ticks=gene_number,
	           labels=gene_ids,
	           rotation=90,
	           fontsize=font_size)
	plt.yticks(fontsize=16)
	plt.ylabel('Samples covered by panel (%)', fontsize=16)
	plt.ylim(50, 100)
	plt.tight_layout()
	plt.savefig(args.output + output_modifier + 'Set_coverage_' + str(datetime.datetime.now().strftime('%Y-%m-%d')) + '.pdf', format='pdf')
	plt.close()


def clusterplot(args, data_matrix, output_modifier=''):
	"""
    Creates clustered mutation matrix based on target genes.
    'data_matrix' is a binary array as a pd.DataFrame; columns=sample, rows=genes. '1' indicates that
    gene X is mutated in sample Y, '0' indicates that it is not.
    """
	if args.numGenes > 50:
		args.numGene = 50
	if len(data_matrix.index) > args.numGenes:
		data_matrix = data_matrix.head(args.numGenes)
	
	if args.targets and not args.greedy_coverage:
		data_matrix = data_matrix.sort_index()
	
	data_matrix.to_csv(args.output + output_modifier + 'Mutation_matrix_datafile_' + str(datetime.datetime.now().strftime('%Y-%m-%d')) + '.tsv', sep='\t')
	
	if args.fixed_matrix_width:
		column_width = 0.025
		height = len(data_matrix) / 2
		plt.figure(figsize=(column_width * data_matrix.shape[1], height))
	else:
		height = len(data_matrix) / 2
		plt.figure(figsize=(8, height))
	
	sns.clustermap(data_matrix,
	               row_cluster=False,
	            cmap='plasma',
	            cbar=False)
	plt.xticks([])
	plt.yticks(rotation=0, fontsize=16)
	plt.tight_layout()
	
	plt.savefig(args.output + output_modifier + 'Mutation_matrix_' + str(datetime.datetime.now().strftime('%Y-%m-%d')) + '.pdf', format='pdf')
	plt.close()


def joint_occurrence_matrix(args, data_matrix, output_modifier=''):
	"""
   Creates a joint-occurrence matrix between genes expressed as a percentage. Allows for a visual identification of which genes co-occur.
   'data_matrix' is a binary array as a pd.DataFrame; columns=sample, rows=genes. '1' indicates that
   gene X is mutated in sample Y, '0' indicates that it is not.
   """
	if args.numGenes > 50:
		args.numGene = 50
	if len(data_matrix.index) > args.numGenes:
		data_matrix = data_matrix.head(args.numGenes)
	
	draw_annotations = True
	if data_matrix.shape[0] >= 20:
		draw_annotations = False
		jo_value_size = 0
	elif data_matrix.shape[0] > 10:
		jo_value_size = 6
	else:
		jo_value_size = 12
	
	max_font_size = 14
	font_size = min(max_font_size, 200 / len(data_matrix))
	
	ind = list(data_matrix.index)
	co_occurrence_matrix = data_matrix.dot(data_matrix.T)
	co_occurrence_matrix_diagonal = np.diagonal(co_occurrence_matrix)
	with np.errstate(divide='ignore', invalid='ignore'):
		co_occurrence_matrix_percentage = np.nan_to_num(np.true_divide(co_occurrence_matrix, co_occurrence_matrix_diagonal))
	df = pd.DataFrame(data=co_occurrence_matrix_percentage, columns=ind, index=ind)
	df.to_csv(args.output + output_modifier + 'Joint_occurence_matrix_datafile_' + str(datetime.datetime.now().strftime('%Y-%m-%d')) + '.tsv', sep='\t')
	plt.figure(figsize=(8, 6))
	sns.heatmap(df,
	            cmap='viridis',
	            xticklabels=True,
	            yticklabels=True,
	            annot=draw_annotations,
	            annot_kws={'size': jo_value_size},
	            vmin=0,
	            vmax=1)
	plt.yticks(rotation=45, fontsize=font_size)
	plt.xticks(rotation=45, fontsize=font_size)
	plt.savefig(args.output + output_modifier + 'Joint_occurrence_matrix_' + str(datetime.datetime.now().strftime('%Y-%m-%d')) + '.pdf', format='pdf')
	plt.close()


def plot_column_sums(args, data_matrix, output_modifier=''):
	"""
    Creates a histogram of co-occurring mutation gene occurrences from only genes selected through either a targets list or top genes from all samples.
    'data_matrix' is a binary array as a pd.DataFrame; columns=sample, rows=genes. '1' indicates that
    gene X is mutated in sample Y, '0' indicates that it is not.
    """
	if args.numGenes > 50:
		args.numGene = 50
	if len(data_matrix.index) > args.numGenes:
		data_matrix = data_matrix.head(args.numGenes)
	
	totals = data_matrix.sum()
	
	plt.figure(figsize=(8, 6))
	sns.histplot(x=totals,
	             discrete=True,
	             color='green',
	             edgecolor='black',
	             legend=False,
	             alpha=0.2)
	plt.ylabel('Frequency', fontsize=16)
	plt.xlabel('Target variants in sample', fontsize=16)
	plt.locator_params(axis="both", integer=True, tight=True)
	plt.yticks(fontsize=16)
	plt.xticks(fontsize=16)
	plt.savefig(args.output + output_modifier + 'Filtered_genes_mutations_per_sample_' + str(datetime.datetime.now().strftime('%Y-%m-%d')) + '.pdf', format='pdf')
	plt.close()


def total_sample_mutations(args, maf_list, output_modifier=''):
	"""
    Creates a histogram of total non-synonymous gene mutations per sample.
    'maf_list' is the list of MAF files to be examined. MAF files should be in the working directory.
    Assumes all samples have at least one mutation.
    """
	
	counts = []
	zero_variant_samples = 0
	for file in maf_list:
		with open(file, 'rt') as f:
			total_mutations = [line.split('\t')[0] for line in f if not line.startswith('#')]
		counts.append(len(total_mutations))
	for item in counts:
		if item == 0:
			zero_variant_samples += 1
	counts = [x for x in counts if x != 0]

	binwidth = 0.0625
	bins = np.arange(np.floor(0), np.ceil(6) + binwidth, binwidth)
	
	plt.figure(figsize=(8, 6))
	sns.histplot(x=counts,
	             log_scale=True,
	             color='green',
	             edgecolor='black',
	             alpha=0.2,
	             bins=bins)
	plt.ylabel('Frequency', fontsize=16)
	plt.xlabel('Number of Variants', fontsize=16)
	plt.yticks(fontsize=16)
	plt.xticks(fontsize=16)
	plt.legend([f'Zero variant samples: {zero_variant_samples}'],
	           loc='upper right',
	           handlelength=0,
	           handleheight=0, fontsize=14)
	plt.xlim(10e0, 10e5)
	plt.savefig(args.output + output_modifier + 'Total_mutations_per_sample_' + str(datetime.datetime.now().strftime('%Y-%m-%d')) + '.pdf', format='pdf')
	plt.close()


def greedy_set_cover(data_matrix, size):
	"""
    Greedy algorithm for the set coverage problem. Takes the most common gene across each sample, followed by
    selection of the gene which is present in the most samples not covered by the original. The process is
    continued until complete set coverage.

    'data_matrix' is a binary array as a pd.DataFrame; columns=sample, rows=genes. '1' indicates that
    gene X is mutated in sample Y, '0' indicates that it is not.
    'size' is an integer and is used as the maximum number of genes to iterate through during greedy set coverage
    """
	
	# if len(data_matrix.index) > 100:
	#     data_matrix = data_matrix.head(100)
	
	cover = pd.DataFrame(data_matrix.iloc[[0]])
	s = [i for i, c in enumerate(cover.columns) if 0.0 in cover[c].values]
	df2 = data_matrix[s]
	for _ in range(1, size + 1):
		df2 = df2[s]
		df2['Count'] = df2.sum(axis=1)
		df2 = df2.sort_values(by=['Count'], ascending=False)
		df2 = df2.iloc[:, :-1]
		index = df2.index[0]
		cover.loc[index] = data_matrix.loc[index]
		column_sum = cover.sum().to_list()
		s = [ind for ind, x in enumerate(column_sum) if x == 0.0]
	return cover


def hierarchial_clustering(args, data_matrix, sample_mapping, num):
	"""
    Creates an ordered hierarchical cluster between samples. Allows for a visual identification similarity between samples.
    Will also separate out clusters for individual cluster analysis: provide and INT value for the number of branches desired.
    'data_matrix' is a binary array as a pd.DataFrame; columns=sample, rows=genes. '1' indicates that
    gene X is mutated in sample Y, '0' indicates that it is not.
    """
	if args.numClusterGenes < args.numGenes:
		args.numClusterGenes = args.numGenes
	if len(data_matrix.index) > args.numClusterGenes:
		data_matrix = data_matrix.head(args.numClusterGenes)
	
	data_matrix = data_matrix.T
	linkage_data = linkage(data_matrix, method='ward', metric='euclidean', optimal_ordering=True)
	dendrogram(linkage_data)
	plt.savefig(args.output + 'Dendrogram_' + str(datetime.datetime.now().strftime('%Y-%m-%d')) + '.pdf', format='pdf')
	plt.close()
	
	dct = dict([(i, {i}) for i in range(data_matrix.shape[0])])
	container = [dct]
	for i, row in enumerate(linkage_data, data_matrix.shape[0]):
		dct[i] = dct[row[0]].union(dct[row[1]])
		del dct[row[0]]
		del dct[row[1]]
		container.append(list(dct.values()))
	
	sub_clusters = container[-num]
	
	clusters_dict = {}
	for i, branch in enumerate(sub_clusters, start=1):
		cluster_sample_ids = []
		for column_number in branch:
			if column_number in sample_mapping:
				cluster_sample_ids.append(sample_mapping[column_number])
		cluster_number = 'Cluster_' + str(i)
		clusters_dict[cluster_number] = cluster_sample_ids
		os.mkdir(args.output + cluster_number)
		
		branch = list(branch)
		sub_df = data_matrix.T[branch]
		
		sub_df['Count'] = sub_df.sum(axis=1)
		sub_df = sub_df.sort_values(by=['Count'], ascending=False)
		sub_df = sub_df.iloc[:, :-1]
		
		clusterplot(args, sub_df, cluster_number + '/')
		gene_set_coverage(args, sub_df, cluster_number + '/')
		joint_occurrence_matrix(args, sub_df, cluster_number + '/')
		plot_column_sums(args, sub_df, cluster_number + '/')
		total_sample_mutations(args, cluster_sample_ids, cluster_number + '/')
	
	max_length = max(len(lst) for lst in clusters_dict.values())
	for key, value in clusters_dict.items():
		if len(value) < max_length:
			clusters_dict[key] += [None] * (max_length - len(value))
	cluster_samples_tsv = pd.DataFrame(clusters_dict)
	cluster_samples_tsv.to_csv(args.output + 'Clustered_samples.tsv', sep='\t', index=False)
