"""
@date: 28rd May 2022
@author: Anna Ketteler
"""

import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), 'GENIE3'))
sys.path.append(os.path.join(os.path.dirname(__file__)))
from GENIE3_python import GENIE3 as gn
import preprocessing as prp
import argparse
import numpy as np
import os
import pickle
import pandas as pd
import csv
import datetime
from sklearn.preprocessing import Normalizer

# partitioning of data set based on confounder
# required: preprocessed data; one col per gene, one per sample; sample array in same order; gene array in same order
# specify data set path, genes you want to take into account and confounder for partitioning

parser = argparse.ArgumentParser(prog='data_partitioning.py',
                                description='CLI for running GENIE3 on confounder-based vs random partitions.',
                                usage='%(prog)s [options] path',)

parser.add_argument('cancer_type',
                       metavar='cancer_type',
                       type=str,
                       help='TCGA cancer type abbreviation (e.g. ´TCGA-BLCA´)')
parser.add_argument('confounder',
                       metavar='confounder',
                       type=str,
                       help='Confounder: ´sex´, ´race´, ´ethnicity´')
parser.add_argument('n',
                       metavar='n',
                       type=int,
                       help='Number of random partitions to generate')
parser.add_argument('m',
                       metavar='m',
                       type=int,
                       help='Number of confounder based partitions to generate')
parser.add_argument('nthreads',
                       metavar='nthreads',
                       type=int,
                       help='Number of threads used by GENIE3')
parser.add_argument('ntrees',
                       metavar='ntrees',
                       type=int,
                       help='Number of trees in ensemble trees generated by GENIE3')
parser.add_argument('-prp',
                        dest='do_preprocessing',
                        action='store_true', # stores the Boolean value True when the corresponding optional argument is specified
                        help='Preprocess the raw data. This option is required at the very first start of the script.')


####################################################################
# parameters
####################################################################
args = parser.parse_args()

cancer_type = args.cancer_type
confounder = args.confounder
do_preprocessing = args.do_preprocessing
n = args.n
m = args.m
nthreads=args.nthreads
ntrees=args.ntrees

####################################################################
# path handling
####################################################################
dateStart = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
print('Result folder time stamp: ' + str(dateStart))
cwd = os.getcwd()
raw_data_dir = os.path.join(cwd, 'input', str(cancer_type)+'_raw_data')
data_dir = os.path.join(cwd, 'input', str(cancer_type)+'_data')
results_dir = os.path.join(cwd, 'results')
timestamp_dir = os.path.join(results_dir, str(dateStart))
conf_results_dir = os.path.join(timestamp_dir, str(cancer_type)+'_conf_results')
rnd_results_dir = os.path.join(timestamp_dir, str(cancer_type)+'_rnd_results')

if not os.path.exists(results_dir):
    os.mkdir(results_dir)
if not os.path.exists(timestamp_dir):
    os.mkdir(timestamp_dir)
if not os.path.exists(conf_results_dir):
    os.mkdir(conf_results_dir)
if not os.path.exists(rnd_results_dir):
    os.mkdir(rnd_results_dir)

####################################################################
# start of script
####################################################################

if do_preprocessing:

    ######################
    # load data, cut column with gene names, header with sample names and then transpose the data matrix
    ######################
    data = np.genfromtxt(fname=os.path.join(raw_data_dir, str(cancer_type)+'.htseq_fpkm.tsv'), delimiter="\t", skip_header=1) # cut header with sample identifiers
    data = data[:,1:]
    (rows, cols) = data.shape
    data = np.transpose(data) # transpose: now one col per gene, one row per sample

    # read gene names and sample ids in separate arrays
    gene_identifiers = np.empty(rows, dtype=object)
    sample_identifiers = np.empty(cols, dtype=object)
    i = -1
    with open(os.path.join(raw_data_dir, str(cancer_type)+'.htseq_fpkm.tsv'),encoding='utf8') as tsvfile:
        tsvreader = csv.reader(tsvfile, delimiter="\t")
        for line in tsvreader:
            if i == -1:
                sample_identifiers = np.array(line[1:], dtype=object)
            if i > -1:
                gene_identifiers[i] = str(line[0])
            i = i + 1

    if not os.path.exists(data_dir):
        os.mkdir(data_dir)

    # save gene identifiers, sample identifiers and data
    pickle.dump(sample_identifiers, open(os.path.join(data_dir, str(cancer_type)+'_sample_identifiers.pkl'), 'wb'))
    pickle.dump(gene_identifiers, open(os.path.join(data_dir, str(cancer_type)+'_gene_identifiers.pkl'), 'wb'))
    pickle.dump(data, open(os.path.join(data_dir, str(cancer_type)+'_data_set.pkl'), 'wb'))
    
    ######################
    # create dictionary for translation of gene identifiers to gene names with mapping file provided by tcga
    ######################
    """
    gene_dict_path = os.path.join(cwd, 'TCGA-BLCA_raw_data', 'gencode.v22.annotation.gene.probeMap')
    gene_dict_data = np.genfromtxt(fname=gene_dict_path, delimiter="\t", skip_header=1, dtype=str)

    ident_list = list(gene_dict_data[:, 0])
    gene_name_list = list(gene_dict_data[:, 1])
    zipped = zip(ident_list, gene_name_list)
    gene_dict = dict(zipped)

    path = os.path.join(cwd, 'TCGA-BLCA_data', 'TCGA-BLCA_gene_dict.pkl')
    with open(path, 'wb') as handle:
        pickle.dump(gene_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
    """

####################################################################
# prepare final data set
####################################################################

######################
# get data: one col per gene, one row per sample
######################
with open(os.path.join(data_dir, str(cancer_type)+'_data_set.pkl'), 'rb') as f:
    raw_data = pickle.load(f)

######################
# get gene identifiers and remove the version identifier and the dot
######################
with open(os.path.join(data_dir, str(cancer_type)+'_gene_identifiers.pkl'), 'rb') as f:
    all_genes = pickle.load(f)
all_genes = np.array([gene.split('.')[0] for gene in all_genes])

######################
# get sample identifiers
######################
with open(os.path.join(data_dir, str(cancer_type)+'_sample_identifiers.pkl'), 'rb') as f:
    all_samples = pickle.load(f)

######################
# recombine into pandas data frame
######################
df = pd.DataFrame(raw_data, columns=all_genes, index=all_samples)
gene_names = np.array(all_genes)
for col in df:
    if df[col].sum() == 0:
            df = df.drop(col, axis=1)
            gene_names = np.delete(gene_names, np.where(gene_names == col))
gene_names = np.array(df.columns.copy())

######################
# filter samples (rows): only keep 'Primary Tumor' samples
######################
tum_samples = prp.getPrimaryTumorIndices(os.path.join(raw_data_dir, str(cancer_type)+'.GDC_phenotype.tsv'))[1]
df = df.filter(items = tum_samples, axis=0)

######################
# get list of regulators
######################
regulators = np.genfromtxt(fname=os.path.join(raw_data_dir, 'TFs_Ensembl_v_1.01.txt'), delimiter="\t", dtype=str)

######################
# group sample ids based on confounder-induced partitions
######################
res = prp.getSampleIDsByConfounder(data_dir, os.path.join(raw_data_dir, str(cancer_type)+'.GDC_phenotype.tsv'), confounder) # returns the confounder based class names and an arrays of sample ids corresponding each to one class
classes = res[0]
conf_partition = res[1]

####################################################################
# start GENIE3 on confounder based partitions
####################################################################

for k in range(m):
    i = 0
    for block in conf_partition:

        if len(block) <=1:
            print('Block must contain at least two elements, otherwise, normalization to unit variance will cause NaN values.')
            break

        # generate partition: select samples from data set based on the index set that was prepared before
        final_df = df.filter(items = block, axis='index')

        gene_names_cpy = gene_names.copy()
        final_df, gene_names_cpy = prp.normalizeToUnitVariance(final_df, gene_names_cpy)
        gene_names_cpy = np.array(final_df.columns.copy())

        # filter from the set of regulators such genes that are not present in the filtered data set 
        mask = np.isin(regulators, gene_names_cpy)
        regulators = np.array(regulators)[mask]
        
        # only keep gene columns in the data frame that the user wishes to examine
        final_df = final_df[gene_names_cpy]
        # make sure that the gene_names_cpy list matches the cols of the data frame again
        mask = np.isin(gene_names_cpy, final_df.columns.values)
        gene_names_cpy = np.array(gene_names_cpy)[mask]

        gene_names_cpy = gene_names_cpy.tolist()
        regulators = regulators.tolist() # TODO: print error message if there are no genes in regulators left and skip iteration
        if len(regulators) < 1 or len(gene_names_cpy) < 1:
            print('no overlap of regulators and genes in data set. Skipping iteration.')
            break
        data = final_df.to_numpy()
        
        p = len(gene_names_cpy)

        # Use Random Forest method
        tree_method='RF'
        # Number of randomly chosen candidate regulators at each node of a tree
        K = p-1
        # Number of trees per ensemble
        #ntrees = 50
        # Run the method with these settings
        VIM3 = gn.GENIE3(data,tree_method=tree_method,K=K,nthreads=nthreads,ntrees=ntrees)
        print(VIM3)

        gn.get_link_list(VIM3,gene_names=gene_names_cpy, regulators=regulators,file_name=os.path.join(conf_results_dir, 'ranking_iter_'+str(k)+'_confBlock_'+str(classes[i])+'.txt'))

        i = i+1

####################################################################
# start GENIE3 on n random partitions with corresponding fractions
####################################################################

# compute fractions for blocks in random partitions
pos = 0
cuts = []
for block in conf_partition:
    cuts.append(pos+len(block))
    pos = pos + len(block)

print(cuts)

# generate partitions according to the fractions computed above
for k in range(n): # n random partitions
    partition = []
    final_df = df.copy()
    samples_cpy = (final_df.index.values).copy()
    np.random.shuffle(samples_cpy)
    partition = np.split(samples_cpy, cuts[:-1])
    i = 0
    # iterate through partition
    for block in partition:

        final_df = df[df.index.isin(block)]

        gene_names_cpy = gene_names.copy()
        final_df, gene_names_cpy = prp.normalizeToUnitVariance(final_df, gene_names_cpy)
        gene_names_cpy = np.array(final_df.columns.copy())

        mask = np.isin(regulators, gene_names_cpy)
        regulators = np.array(regulators)[mask]
        
        # only keep gene columns in the data frame that the user wishes to examine
        final_df = final_df[gene_names_cpy]

        # make sure that the gene_names_cpy list matches the cols of the data frame again
        mask = np.isin(gene_names_cpy, final_df.columns.values)
        gene_names_cpy = np.array(gene_names_cpy)[mask]

        gene_names_cpy = gene_names_cpy.tolist()
        regulators = regulators.tolist()
        data = final_df.to_numpy()
        if len(regulators) < 1 or len(gene_names_cpy) < 1:
            print('no overlap of regulators and genes in data set. Skipping iteration.')
            break

        p = len(gene_names_cpy)

        # Use Random Forest method
        tree_method='RF'
        # Number of randomly chosen candidate regulators at each node of a tree
        K = p-1
        # Number of trees per ensemble
        #ntrees = 50
        # Run the method with these settings
        VIM3 = gn.GENIE3(data,tree_method=tree_method,K=K,nthreads=nthreads,ntrees=ntrees)

        gn.get_link_list(VIM3,gene_names=gene_names_cpy, regulators=regulators,file_name=os.path.join(rnd_results_dir, 'ranking_iter_'+str(k)+'_rndBlock_'+str(classes[i])+'.txt'))

        i = i+1

print('Result folder time stamp: ' + str(dateStart))
