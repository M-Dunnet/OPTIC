# OPTIC (Oncogene Panel Tester for Identifying Cancer)
Oncogene Panel Tester for the Identification of Cancer (OPTIC) is a computational pipeline for determining suitable target genes for cancer panel creation.

## Installation
OPTIC can be downloaded using the following command:
`git clone https://github.com/M-Dunnet/OPTIC`

### Dependencies

* contourpy==1.1.0
* cycler==0.11.0
* fonttools==4.42.1
* importlib-resources==6.0.1
* kiwisolver==1.4.5
* matplotlib==3.7.2
* numpy==1.25.2
* packaging==23.1
* pandas==2.1.0
* Pillow==10.0.0
* pyparsing==3.0.9
* python-dateutil==2.8.2
* pytz==2023.3
* scipy==1.11.2
* seaborn==0.12.2
* six==1.16.0
* tzdata==2023.3
* zipp==3.16.2

Dependencies can be installed with the requirements.txt file. To install, first create a virtual environment:

    python3 -m venv venv_OPTIC

To activate the virtual environment use the command navigate to the environment directory and use the command:

    source venv_OPTIC/bin/activate

With the environment active, install the dependencies:

    pip install -r requirements.txt
    

## How OPTIC works

OPTIC aids the generation of small sequencing panels. OPTIC takes MAF files as input and creates a binary array that indicates the presence of somatic mutations in every gene for each sample. The pipeline is designed to operate through a series of iterations, with each step followed by manual inspection. By default, each iteration of OPTIC produces several plots and data files to aid in panel development. OPTIC will produce a mutation matrix plot detailing the gene mutations of each sample, a co-occurrence plot detailing which genes are commonly mutated together, and a set coverage plot. The set coverage plot shows the increase in the number of cumulative samples captured by the addition of individual genes. OPTIC also produces a histogram of the total number of variants per sample. Total gene mutation frequency, gene variant counts (i.e. the number of times each specific variant was found in the dataset) and gene mutation type counts (i.e. the number of missense, nonsense, deletion, and insertion mutations per gene) are also reported in the output.

The OPTIC workflow (Figure 1) begins with preparing a variant filter file to remove non-pathogenetic variants. Next, an iteration of OPTIC can be run with any of the following options: 
1. Default, which will examine all genes from all samples. 
2. Hierarchical clustering, which will use Ward's minimum variance method and Euclidean distance to cluster samples based on mutation profiles.
3. Greedy set coverage, which will choose the fewest number of genes to cover the highest number of samples.
4. Targeted, in which the user will provide a predefined list of genes for OPTIC to examine.

OPTIC was designed such that hierarchical clustering was performed first to separate tumours with different molecular profiles. For example, separating hypermutated and non-hypermutated cancers so that analysis of each can be performed separately. This is followed by the greedy set coverage option to find the most efficient set of genes for a sequencing panel, and then by the targeted option to assess the panel. Nevertheless, it is not a requirement that OPTIC is used in this order. Full descriptions of each option are provided below.

### Generating a filter file
The purpose of the filter is to remove passenger mutations and variants that are not clinically relevant to the disease type. The filter file is produced using data from the Cancer Mutation Census from the Catalogue of Somatic Mutations in Cancer (COSMIC) (note this is not the census gene mutations file found in the 'COSMIC' project, but the data from the 'Cancer Mutation Census' project). The `Generate_filter.py` script transforms the cancer mutation census into a format readable by OPTIC, and only variants within this file are analysed. OPTIC offers the ability to exclude non-pathogenic variants as determined via ClinVar classification and/or to filter based upon COSMIC’s variant classification. The resulting filter can be directly used by OPTIC or can be further modified with the `Hotspot_Genes.py` script to remove variants outside of mutation hotspot locations in specific genes. For example, when targeting the BRAF V600E variant (chr7: 140753336, GRCh38), all BRAF variants outside of this position will be excluded, while all other genes will be unaffected. 

During filtering, OPTIC examines both the genomic position of the mutation, as well as the amino acid change. Only one match in either category is required. This prevents slight discrepancies in mapping algorithms from unintentionally removing relevant mutations. For example, aligners using either 0- or 1-based indexing, or small indels in repeat regions, particularly homopolymer repeats, can create alignment ambiguity. This results in multiple ways to represent the mutation despite the functional change being identical. Likewise, the position of an amino acid change can also be represented in multiple ways, depending on the reference transcript used in variant calling. The transcript ID for variants used by OPTIC is in the filter file, listed in column two. If no filter file is provided, all mutations are considered.

### Hierarchical clustering
When using hierarchical clustering, the file list will be separated into multiple clusters based on molecular subtype; for example, hypermutated and non-hypermutated colorectal cancers. Clustering enables groups of samples with different mutation profiles to be examined separately, allowing for more precise identification of unique patterns within each group. Furthermore, it may enable the detection of low-frequency subtypes driven by specific genes not easily identifiable in large sample sizes. Clusters are generated based on somatic mutations within the binary mutation array using Ward's minimum variance method and Euclidean distance as the measure of dissimilarity. The number of generated clusters and number of genes used in clustering can be set by the user, however, a low number of genes can introduce distance ties, which are broken randomly. Therefore, it is recommended to use as many genes as possible. All plots and data files are generated for each generated cluster. 

### Greedy Set Coverage Algorithm
Finding the minimum number of genes required to cover the maximum number of CRC samples is an example of the set cover problem. Here, the universal set U represents CRC samples, and each subset S1, S2,…,Sn corresponds to a gene that is mutated in specific samples within U. The objective is to find the smallest number of genes (subsets) whose mutations collectively cover all the CRC samples in U. This approach minimizes panel size and sequencing requirements while maintaining maximum sample coverage. OPTIC uses a greedy algorithm to address this problem. It begins by selecting the gene with the highest mutation rate across all samples. Subsequent genes are chosen based on their ability to cover the maximum number of previously uncovered samples. 

When the option for set coverage is used, the order of genes in the mutation matrix and set coverage plot is determined by the number of new samples covered, rather than by total mutation frequency. In some cases, certain samples have unique sets of variants not shared with any other samples. If these samples contain multiple mutations, any one of these mutations could represent the sample. However, the specific gene used for coverage is selected randomly. When multiple samples have this characteristic, they appear as a long tail of single genes covering individual samples in the set coverage plots. It's important to note that these genes are sample-specific, their order is chosen randomly, and they may not be the only genes capable of covering the sample. 

### Targeted OPTIC
Users can provide OPTIC with a text file containing a list of genes using HUGO gene nomenclature to limit the analysis to that gene set. Targeted OPTIC can be run in conjunction with the greedy set coverage algorithm. Using both options enables the evaluation of initial panels and can show the effect of adding or removing specific genes. For instance, if a gene of interest was not selected by the set cover algorithm, users can manually add it to the panel to assess its influence on the analysis results. Conversely, if a gene included by the set coverage algorithm is found to be undesirable or irrelevant, users have the flexibility to remove it from the panel. This iterative process allows users to refine their gene panels, ensuring that both desired and undesired genes are optimally managed to achieve the most accurate and relevant gene set.

### Additional Scripts. 
Additional scripts are also included to supplement the OPTIC workflow. 
1. The `move_files.sh` script can move files from different subsets (as determined by hierarchical clustering) into separate directories.
2. If the data has been downloaded from c.BioPortal as a aggregated meta-file, The `Extract_from_cbioportal.py` separate each sample into its own file.

### OPTIC directory tree

```
OPTIC
│   README.md
│   OPTIC.py
|   OPTIC_config.txt
|   requirements.txt
│
└───lib
│   │   Anaylsis_functions.py
|   |   Data_functions.py
|   |   Config.py
│   
└───bin
|   │   Generate_filter.py
|   │   Hotspot_genes.py
|   |   move_files.sh
|   |   Extract_from_cbioportal.py
|
└───files
|   |   Sample_targets.txt
|   |   Example_Hotspot_regions.txt
```

## OPTIC.py

OPTIC.py takes MAF files as input. The structure of the MAF is available in the GDC MAF specification along with descriptions for each field. OPTIC will generate a binary mutation array comprising all genes from all samples for analysis. 

If the exact structure of the input MAF files does not exactly match that of the GDC MAF specification, it can be edited within the config file. The config file specifies which columns OPTIC will look at for each set of information. Exact instructions are explained within the config file

### Command-line Arguments
#### Required Arguments

    --input, -i    # The directory location containing input maf files
    --output, -o   # The directory location where files will end up. If nothing provided, output files will be placed in the input folder

#### Analysis Options
    --targets, -t            # Path string to the location of the target file. Only loci in the target file will be analysed
    --greedy_coverage        # Boolean option to run the greedy set coverage algorythm. This will select the subsets of a given gene to cover as many samples as possible, prioritizing subsets that cover the most uncovered elements at each step.
    --hierarchial_clustering # Sets OPTIC to cluster samples by mutation profile. The integer provided sets the number of clusters to analyse

#### Filter Options
    --cosmic_mutants       # Path string to the location of the COSMIC filter file. Only mutations within the filter file will be analysed.
    --numGenes, -n         # Maximum number of genes to use in the analysis. Defaults to 25, maximum of 50`
    --numClusterGenes, -c  # Maximum number of genes to examine during hierarchical clustering. Defaults to 250. Note this number cannot be lower than numGenes.`
    --save_zeros, -z       # Boolean option to save a list of all samples that do not contain mutations in the gene set`
    --use_cds              # Boolean option, if set OPTIC will use CDS changes rather than amino acid changes`

#### Plotting options
    --fixed_matrix_width, -f   # Boolean option to plot the mutation matrix where each sample has a fixed width. If not set, column widths will be scaled relative to the sample number. Default is False
    --plot_total_variants      # Boolean option to plot the total number of variants per sample (includes passenger mutations). Default is False

### Example Usage:
```
python3 OPTIC.py -i <input_location> -o <output_location> --greedy_coverage --cosmic_mutants <filter_file_location> --plot_total_variants
```

### Output files
#### Dendrogram
This is only produced when the hierarchical clustering option is used. It shows the Euclidean distance between each sample based on somatic mutation profiles (including passenger mutations) as determined by Ward's minimum variance. 


#### Mutation matrix
The mutation matrix provides a visual identification of which genes are mutated at specific frequencies. Samples are represented by columns, and genes across rows. Yellow indicates at least one mutation in the specified gene, while blue indicates a wild-type genotype.


#### Cumulative coverage plots
The cumulative coverage plot shows the proportion of samples encapsulated by the addition of different genes. Genes are ordered by total mutation frequency, or, if the set coverage algorithm is used, by total sample coverage increase.

#### Joint occurrence matrix
A joint-occurrence matrix between genes expressed as a percentage. Allows for a visual identification of which genes co-occur. Each cell at the intersection of a row (Gene A) and a column (Gene B) shows the co-occurrence percentage, representing the proportion of samples with a mutation Gene A that also contain a mutation in Gene B.

#### Filtered gene mutations per sample
Histogram of the number of variants per sample after variant filtering. This only includes genes from either the `--targets` list (if used) or the top number of genes equal to the `--numGenes` argument

#### Total sample mutations
A histogram of the total number of mutations per sample. This includes all mutations, regardless if they are in the filter file or not. 

#### Gene mutation frequencies
Tab-delimited list showing the gene name (HUGO nomenclature) and the overall mutation frequency across all samples. 

#### Variant counts file
Tab-delimited list showing the gene name (HUGO nomenclature), specific variant (e.g. c.G12A), and the total number of times it has occurred across all samples. 

#### Variant types file
Tab-delimited list showing the gene name (HUGO nomenclature), the type of mutation (e.g. missense), and the total number of times it has occurred across all samples. 


## Generate_filter.py
Generates a filter file for OPTIC using the COSMIC Cancer Mutation Census dataset. Only variants within the filter file will be considered by OPTIC.

### Arguments

    --input, -i              # File path of the cancer mutation census tsv file
    --output, -o             # Location to save parsed DataFrame. If nothing provided, output files will be placed in the input folder
    --genome_build, -g       # Defines which genome version to use, must be either `GRCh38` (or `hg38`) or `GRCh37` (or `hg19`)
    --filter_t1, -t          # Boolean. If set, all genes not within the tier 1 COSMIC gene census will be excluded. Defaults to False
    --filter_pathogenic, -p  # Boolean. If set, only pathogenic ClinVar mutations will be kept. Defaults to False 
    --filter_cosmic_sig, -c  # Set as INT between 1-4. If set, only mutations at that tier level and higher will be kept (e.g. if set as 2, tier 1 and 2 mutations will be kept. The `Other` tier is renamed to 4. Default is 4 (No filtering)

### Example Usage
```
python3 Generate_filter.py -i <input_location> -o <output_location> -g hg38 -c 3
```

## Hotspot_genes.py
Limits the mutation filter file to only contain variants from specific regions for specified genes.

    --input, -i   # The location of the Parsed COSMIC mutants tsv file
    --output, -o  # Location to save hotspot filtered DataFrame
    --hotspots    # The location of the hotspots file. Must be in TSV format

The hotspots file must contain the gene name in HUGO nomenclature and the genomic start and end positions of the regions of interest. Values must be tab-separated. 

### Example Usage
```
python3 Hotspot_genes.py -i <filter_file_location> -o <output_location> -hotspots <location_of_hotspots_file>
```

## Additional scripts.py
### move_files.sh
Running OPTIC with the hierarchical_clustering will produce a TSV format file showing which samples are in each cluster. `move_files.sh` is a short bash script that can take this file as input, and will move files from an original directory into cluster-specific directories. For example, if the original files are in a directory:

```
Dir containing MAF files
│   File1.maf
│   File2.maf
|   File3.maf
|   File4.maf
```
and the `Clustered_samples.tsv` file produced by OPTIC shows samples were clustered into two groups:
| Cluster 1  | Cluster 2 |
|------------|-----------|
| File1.maf  | File2.maf |
| File3.maf  | File4.maf | 

The `move_files.sh` will create the following directory tree:

```
Output Dir
└───Cluster 1 samples
│   │   File1.maf
|   |   File3.maf
│   
└───Cluster 2 samples
│   │   File2.maf
|   |   File4.maf
```

#### Usage
```
./move_files.sh <tsv_file> <input_directory> <output_directory>"
```
All arguments are required, positional arguments.
* tsv_file: The clustered_samples.tsv file produced by OPTIC
* input_directory: Location of the original directory containing all MAF files
* output_directory: Location where the main output_dir (see above) will be produced.

### Extract_from_cBioPortal.py
 If the data has been downloaded from c.BioPortal as a aggregated meta-file, The `Extract_from_cbioportal.py` separate each sample into its own file using the `Tumor_Sample_Barcode` column.
 
 #### Usage
 ```
 python3 Extract_from_cbioportal.py <filename> <output_location>
 ```
 
    filename:        Location of the meta-file to process
    output_location: Location where split files will end up

Both arguments are positional.
