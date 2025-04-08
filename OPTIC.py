#!/usr/bin/env python3
"""
Oncogene Panel Tester for the Identification of Cancer (OPTIC) is a computational pipeline for determining suitable target genes for cancer panel creation.
"""

import os
import sys
import argparse
import datetime
import warnings

import pandas as pd

PATH = '/'.join(os.path.realpath(__file__).split('/')[:-1]) + '/lib/'
sys.path.append(os.path.abspath(PATH))

from Config import OpticConfigFile
import Data_functions as optic_data
import Analysis_functions as optic_analysis


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
                             'to cover as many samples as possible, prioritizing subsets that cover the most uncovered elements at each step.')
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
                        help='If set, OPTIC will use CDS changes rather than acid changes')
    
    # Plotting Options
    parser.add_argument('--fixed_matrix_width', '-f', default=False, action='store_true',
                        help='Boolean option to plot the mutation matrix where each sample has a fixed width. '
                             'If not set, column widths will be scaled relative to the sample number. Default is False')
    parser.add_argument('--plot_total_variants', default=False, action='store_true',
                        help='Boolean option to plot the total number of variants per sample (includes passenger mutations). Default is False')
    
    return parser.parse_args()


def get_file_info(file):
    raw_sample_variants_dict = {}
    try:
        with open(file, 'rt') as f:
            for j, line in enumerate(f):
                if line.startswith('#'):
                    continue
                
                columns = line.strip().split('\t')

                # Apply column filters if present
                if ('Filter Column 1' in input_file_attributes and
                        'Filter Column 1 Exclusions' in input_file_attributes and
                        columns[input_file_attributes['Filter Column 1']] in input_file_attributes['Filter Column 1 Exclusions']):
                    continue
                
                if ('Filter Column 2' in input_file_attributes and
                        'Filter Column 2 Exclusions' in input_file_attributes and
                        columns[input_file_attributes['Filter Column 2']] in input_file_attributes['Filter Column 2 Exclusions']):
                    continue
                
                variant_key = 'Variant CDS' if args.use_cds else 'Variant Amino Acid'
                raw_sample_variants_dict[j] = [
                    columns[input_file_attributes['Gene']],
                    columns[input_file_attributes['Start_position']],
                    columns[input_file_attributes['Reference Allele']],
                    columns[input_file_attributes['Alternate Allele']],
                    columns[input_file_attributes[variant_key]],
                    columns[input_file_attributes['Variant Type']]
                ]
    except Exception as e:
        print(f"Error processing file {file}: {e}")
    
    return raw_sample_variants_dict


def main(arg, input_file_attributes):
    # Initialize the list of files
    maf_list = []
    for maf_file in os.listdir(args.input):
        if maf_file.endswith(input_file_attributes['File Extension']):
            maf_list.append(maf_file)
    maf_list = list(set(maf_list))
    maf_list.sort()
    
    if len(maf_list) == 0:
        print(f'Error: There are no {input_file_attributes["File Extension"]} files in input folder... '
              f'Check the correct file extension is provided in the config file')
        sys.exit(1)
    
    # Initialize list of target genes if specified
    if args.targets:
        with open(args.targets, 'rt') as f:
            targets = f.read().splitlines()
    
    # Initialize COSMIC mutants filter file if specified
    if arg.cosmic_mutants:
        cosmic_mutants = pd.read_csv(arg.cosmic_mutants, sep='\t')
        cosmic_mutants = cosmic_mutants.astype({'Position': float, 'Variant_CDS': str, 'Variant_AA': str})
        if arg.use_cds:
            cosmic_mutants.rename(columns={'Variant_CDS': 'Variant'}, inplace=True)
        else:
            cosmic_mutants.rename(columns={'Variant_AA': 'Variant'}, inplace=True)
        
        cosmic_mutants['GENOMIC_WT_ALLELE'] = cosmic_mutants['GENOMIC_WT_ALLELE'].fillna('-')
        cosmic_mutants['GENOMIC_MUT_ALLELE'] = cosmic_mutants['GENOMIC_MUT_ALLELE'].fillna('-')
    else:
        cosmic_mutants = None
    
    # Load in file information
    mutation_dict, number = {}, 0
    for file in maf_list:
        number += 1
        raw_sample_variants_dict = get_file_info(file)
        sample_variants_df = pd.DataFrame.from_dict(raw_sample_variants_dict, orient='index')
        sample_variants_df.columns = ['Gene_ID', 'Position', 'Reference_allele', 'Alternate_allele', 'Variant', 'Variant_type']

        gene_dict = optic_data.process_optic_dictionary(args, file, sample_variants_df, cosmic_mutants)
        mutation_dict.update(gene_dict)
        
        print(f"Finished file {number} of {len(maf_list)} ({round(number / len(maf_list) * 100, 2)}%)", ' ' * 20, end='\r')
    
    # Process Variants
    optic_data.count_variant_types(args, mutation_dict)
    if args.targets:
        optic_array, sample_mapping = optic_data.generate_optic_array(mutation_dict, targets)
    else:
        targets = optic_data.get_gene_ids(mutation_dict)
        optic_array, sample_mapping = optic_data.generate_optic_array(mutation_dict, targets)
    
    optic_data.get_total_frequencies(args, mutation_dict, optic_array, maf_list)
    optic_array = optic_array.iloc[:, :-1]
    optic_data.check_optic_array(optic_array)
    
    # Analyse Data
    gsc_optic_array = optic_analysis.greedy_set_cover(optic_array, arg.numGenes)
    
    if args.hierarchial_clustering:
        optic_analysis.hierarchial_clustering(arg, optic_array, sample_mapping, args.hierarchial_clustering)
    
    elif args.greedy_coverage:
        optic_analysis.clusterplot(arg, gsc_optic_array)
        optic_analysis.gene_set_coverage(arg, gsc_optic_array)
        optic_analysis.joint_occurrence_matrix(arg, gsc_optic_array)
        optic_analysis.plot_column_sums(arg, gsc_optic_array)
    
    else:
        optic_analysis.clusterplot(arg, optic_array)
        optic_analysis.gene_set_coverage(arg, optic_array)
        optic_analysis.joint_occurrence_matrix(arg, optic_array)
        optic_analysis.plot_column_sums(arg, optic_array)
    
    if args.plot_total_variants:
        optic_analysis.total_sample_mutations(arg, maf_list)


if __name__ == '__main__':
    warnings.filterwarnings("ignore")
    
    # Define arguments
    args = parse_args()
    
    # Load Configuration file
    config_file = 'OPTIC_config.txt'
    config_reader = OpticConfigFile(config_file)
    config_reader.read_config()
    config_reader.validate_config(args)
    input_file_attributes = config_reader.get_attributes()
    
    # Check arguments are valid
    if not args.input:
        print('Error: an input folder (--input/-i) is required', file=sys.stderr)
        sys.exit(1)
    if not os.path.isdir(args.input):
        print('Error: input (--input/-i) is not a directory', file=sys.stderr)
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
    
    # Make Output locations
    os.chdir(args.input)
    date = str(datetime.datetime.now().strftime('%Y-%m-%d'))
    if args.targets:
        args.output = f'{args.output}Targeted_OPTIC_{date}/'
        if not os.path.exists(args.output):
            os.mkdir(args.output)
    elif args.hierarchial_clustering:
        args.output = f'{args.output}Hierarchial_clustering_{args.hierarchial_clustering}_clusters_{date}/'
        if not os.path.exists(args.output):
            os.mkdir(args.output)
    else:
        args.output = f'{args.output}All_genes_OPTIC_{date}/'
        if not os.path.exists(args.output):
            os.mkdir(args.output)
    
    main(args, input_file_attributes)
