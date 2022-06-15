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
    """Returns the sample ids that contain information obtained from Primary Tumor tissue.
        
        Parameters
        ----------
        pheno_data_path : str
            Path to file which contains phenotypic information. One row per sample, one column per attribute.

        Returns
        -------
        indices : np.array
            Sample IDs of samples from Primary Tumor tissue.
    """
    pheno_data = np.genfromtxt(fname=pheno_data_path, delimiter="\t", dtype=str)
    
    # remove normal and metastatic tissue samples from set of sample ids
    tissue_col_index = np.where(pheno_data[0, :] == 'sample_type.samples')[0][0]
    tum_data = pheno_data[pheno_data[0:, tissue_col_index] == 'Primary Tumor']
    tum_samples = tum_data[:, 0]
    indices = (['primary_tum_samples'], (tum_samples))
    return indices


def getSampleIDsByConfounder(pheno_data_path, confounder):
    """Returns a tuple of two np.arrays with the first containing block string-identifiers of the requested @confounder 
        and the second containing the sample ids corresponding to the blocks.
            
        Parameters
        ----------
        pheno_data_path : str
            Path to file which contains phenotypic information. One row per sample, one column per attribute.

        confounder : str
            Confounder attribute that is to be used to build the partition.
        
        Returns
        -------
        blocks : np.array
            At index i, contains the str-identifier of the block at index i in the @conf_partition.
        
        conf_partition : np.array
            Contains the blocks belonging to the confounder-based partition.
    """
    pheno_data = np.genfromtxt(fname=pheno_data_path, delimiter="\t", dtype=str)
    indices = None

    # split set of sample ids based on confounder expression
    if confounder == 'sex':
        gender_col_index = np.where(pheno_data[0, :] == 'gender.demographic')[0][0]
        female_data = pheno_data[pheno_data[:, gender_col_index] == 'female']
        male_data = pheno_data[pheno_data[:, gender_col_index] == 'male']

        female_samples = female_data[:, 0]
        male_samples = male_data[:, 0]
        (blocks, conf_partition) = (['female', 'male'], [female_samples, male_samples])

    if confounder == 'ethnicity':
        ethn_col_index = np.where(pheno_data[0, :] == 'ethnicity.demographic')[0][0]
        hisp_lat_data = pheno_data[pheno_data[:, ethn_col_index] == 'hispanic or latino']
        nonHisp_nonLat_data = pheno_data[pheno_data[:, ethn_col_index] == 'not hispanic or latino']

        hisp_lat_samples = hisp_lat_data[:, 0]
        nonHisp_nonLat_samples = nonHisp_nonLat_data[:, 0]
        (blocks, conf_partition) = (['hisp_lat', 'non_hisp_lat'], [hisp_lat_samples, nonHisp_nonLat_samples])

    if confounder == 'race':
        race_col_index = np.where(pheno_data[0, :] == 'race.demographic')[0][0]
        asian_data = pheno_data[pheno_data[:, race_col_index] == 'asian']
        african_data = pheno_data[pheno_data[:, race_col_index] == 'black or african american']
        white_data = pheno_data[pheno_data[:, race_col_index] == 'white']

        asian_samples = asian_data[:, 0]
        african_samples = african_data[:, 0]
        white_samples = white_data[:, 0]
        (blocks, conf_partition) = (['asian', 'african', 'white'], [asian_samples, african_samples, white_samples])

    return (blocks, conf_partition)

def normalizeToUnitVariance(df):
    # normalize gene expression data for each gene vector (col) to unit length
    for col in df:
        if df[col].std() != 0:
            df[col] = df[col]/df[col].std()
        else:
            print('df contains column '+str(col) + ' with std()==0. Normalization would produce nan value. Remove column ' + str(col) + ' before normalization.')
    return df

def testDuplicateGenes(genes):
    u, c = np.unique(genes, return_counts=True)
    if len(u)<len(genes):
        print('Cutting version identifier from ensemble ids caused occurrence of duplicates in list of gene identifiers!')

def preprocessToPkl(raw_data_dir, data_dir, cancer_type):

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

def preprocessToCsv(raw_data_dir, data_dir, cancer_type):
    print('start preprocessing...')
    df = (pd.read_csv(os.path.join(raw_data_dir, str(cancer_type)+'.htseq_fpkm.tsv'), sep='\t', header=0, index_col=0)).T
    df = df.rename(columns={'Ensembl_ID':'sample_ID'}) # now the first column contains sample identifiers & should be named accordingly
    if not os.path.exists(data_dir):
        os.mkdir(data_dir)
    print('save pd.DataFrame to csv file ' + str(os.path.join(data_dir, str(cancer_type)+'_data_set.csv')))
    df.to_csv(os.path.join(data_dir, str(cancer_type)+'_data_set.csv'), sep='\t')
