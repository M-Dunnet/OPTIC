import argparse
import pandas as pd
import os


def parse_args():
    parser = argparse.ArgumentParser(description='Separates the files within the meta-file, downloaded from c.bioportal.',
                                     add_help=True,
                                     prefix_chars='-')
    parser.add_argument('filename', type=str, help='The name of the file to process')
    parser.add_argument('output_location', type=str, help='Location where split files will end up')
    return parser.parse_args()


def main(args):
    combined_dataframe = pd.read_csv(args.filename, sep='\t', skiprows=0, low_memory=False)
    dataframe_grouped_by_sample = combined_dataframe.groupby('Tumor_Sample_Barcode')
    for sample_id, group in dataframe_grouped_by_sample:
        group.to_csv(f"{args.output_location}{sample_id}.tsv", sep="\t", index=False)


if __name__ == '__main__':
    args = parse_args()
    if not args.output_location.endswith('/'):
        args.output_location = args.output_location + '/'
    
    if not os.path.exists(args.output_location):
        os.makedirs(args.output_location)
        print(f'Directory {args.output_location} did not exist, and has now been created.')

    main(args)
