import numpy as np
import os
import pickle
import pandas as pd
import csv
from sklearn.preprocessing import Normalizer

####################################################################
# helper methods
####################################################################

def mapEsidToName(esid_list):
    dict_path = os.path.join(cwd, 'TCGA-BLCA_data', 'TCGA-BLCA_gene_dict.pkl')
    with open(dict_path, 'rb') as f:
        gene_dict = pickle.load(f)

    gene_name_list = []
    for i in genes:
        gene_name_list.append(gene_dict.get(i))

    return gene_name_list

def getPrimaryTumorIndices(pheno_data_path):

    pheno_data = np.genfromtxt(fname=pheno_data_path, delimiter="\t", dtype=str)
    
    # remove normal and metastatic tissue samples from set of sample ids
    tissue_col_index = np.where(pheno_data[0, :] == 'sample_type.samples')[0][0]
    tum_data = pheno_data[pheno_data[0:, tissue_col_index] == 'Primary Tumor']

    tum_samples = tum_data[:, 0]
    indices = (['primary_tum_samples'], (tum_samples))

    return indices


def getSampleIDsByConfounder(save_dir, pheno_data_path, confounder):
    pheno_data = np.genfromtxt(fname=pheno_data_path, delimiter="\t", dtype=str)
    
    data_dir = save_dir
    if not os.path.exists(data_dir):
        os.mkdir(data_dir)

    indices = None

    # split set of sample ids based on confounder expression
    if confounder == 'sex':
        gender_col_index = np.where(pheno_data[0, :] == 'gender.demographic')[0][0]
        female_data = pheno_data[pheno_data[:, gender_col_index] == 'female']
        male_data = pheno_data[pheno_data[:, gender_col_index] == 'male']

        pickle.dump(female_data, open(os.path.join(data_dir, 'TCGA-BLCA_female_data.pkl'), 'wb'))
        pickle.dump(male_data, open(os.path.join(data_dir, 'TCGA-BLCA_male_data.pkl'), 'wb'))

        female_samples = female_data[:, 0]
        male_samples = male_data[:, 0]
        indices = (['female', 'male'], (female_samples, male_samples))

    if confounder == 'ethnicity':
        ethn_col_index = np.where(pheno_data[0, :] == 'ethnicity.demographic')[0][0]
        hisp_lat_data = pheno_data[pheno_data[:, ethn_col_index] == 'hispanic or latino']
        nonHisp_nonLat_data = pheno_data[pheno_data[:, ethn_col_index] == 'not hispanic or latino']

        pickle.dump(hisp_lat_data, open(os.path.join(data_dir, 'TCGA-BLCA_hisp_lat_data.pkl'), 'wb'))
        pickle.dump(nonHisp_nonLat_data, open(os.path.join(data_dir, 'TCGA-BLCA_nonHisp_nonLat_data.pkl'), 'wb'))

        hisp_lat_samples = hisp_lat_data[:, 0]
        nonHisp_nonLat_samples = nonHisp_nonLat_data[:, 0]
        indices = (['hisp_lat', 'non_hisp_lat'], (hisp_lat_samples, nonHisp_nonLat_samples))

    if confounder == 'race':
        race_col_index = np.where(pheno_data[0, :] == 'race.demographic')[0][0]
        asian_data = pheno_data[pheno_data[:, race_col_index] == 'asian']
        african_data = pheno_data[pheno_data[:, race_col_index] == 'black or african american']
        white_data = pheno_data[pheno_data[:, race_col_index] == 'white']

        pickle.dump(asian_data, open(os.path.join(data_dir, 'TCGA-BLCA_asian_data.pkl'), 'wb'))
        pickle.dump(african_data, open(os.path.join(data_dir, 'TCGA-BLCA_african_data.pkl'), 'wb'))
        pickle.dump(white_data, open(os.path.join(data_dir, 'TCGA-BLCA_white_data.pkl'), 'wb'))

        asian_samples = asian_data[:, 0]
        african_samples = african_data[:, 0]
        white_samples = white_data[:, 0]
        indices = (['asian', 'african', 'white'], (asian_samples, african_samples, white_samples))

    return indices

def normalizeToUnitVariance(df, gene_names):
    # normalize gene expression data for each gene vector (col) to unit length
    gene_names = np.array(gene_names)
    for col in df:
        if df[col].std() != 0: # this is unnecessary, .std() would be sufficient, but I have to ask how to treat genes without variation
            df[col] = df[col]/df[col].std()
        else:
            df = df.drop(col, axis=1)
            gene_names = np.delete(gene_names, np.where(gene_names == col))
            print('removed gene ' + str(col) + ' from genes that are taken into account, because there were only 0 entries')
    return df, gene_names

def testDuplicateGenes(genes):
    u, c = np.unique(genes, return_counts=True)
    return len(u) < len(genes)

def preprocess(raw_data_dir, cancer_type):

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
