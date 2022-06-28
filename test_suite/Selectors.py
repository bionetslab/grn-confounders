from enum import Enum
from .GENIE3Wrapper import GENIE3Wrapper
from .ARACNEWrapper import ARACNEWrapper
from .LSCONWrapper import LSCONWrapper
import pandas as pd
import numpy as np
import os

class AlgorithmSelector(Enum):
    """Enum specifying which network inference algorithm should be used."""
    ARACNE = 'ARACNE'
    GENIE3 = 'GENIE3'
    #LSCON = 'LSCON'

    def __str__(self):
        return self.value

class CancerTypeSelector(Enum):
    """Enum specifying which cancer type should be investigated."""
    BLCA = 'BLCA'

    def __str__(self):
        return self.value

class ConfounderSelector(Enum):
    """Enum specifying the confounder whose effect is to be examined."""
    SEX = 'sex'
    RACE = 'race'
    AGE = 'age'

    def __str__(self):
        return self.value

def get_algorithm_wrapper(algorithm_selector):
    """Returns the appropriate algorithm based on the selection.

    Parameters
    ----------
    algorithm_selector : AlgorithmSelector
        Specifies which algorithm should be used.
    """
    if algorithm_selector == AlgorithmSelector.GENIE3:
        return GENIE3Wrapper()
    elif algorithm_selector == AlgorithmSelector.ARACNE:
        return ARACNEWrapper()
    #elif algorithm_selector == AlgorithmSelector.LSCON:
        #return LSCONWrapper()

def download_TCGA_expression_data(cancer_type_selector):
    """Saves TCGA gene expression RNAseq - HTSeq - FPKM data for the specifies @cancer_type obtained from UCSC Xena in /data.

    Parameters
    ----------
    cancer_type_selector : CancerTypeSelector
        Specifies the cohort that the phenotype file is to be downloaded for.
    """
    cwd = os.getcwd()
    url = ""
    if cancer_type_selector == CancerTypeSelector.BLCA:
        url = "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-BLCA.htseq_fpkm.tsv.gz"
    df = pd.read_csv(url, delimiter='\t', index_col='Ensembl_ID').T
    #for col in df:
        #df.rename({col: col.split('.')[0]}, axis=1, inplace=True) # TODO
    df.to_csv(os.path.join(cwd, 'data', 'TCGA-'+str(cancer_type_selector)+'.htseq_fpkm.tsv'), sep='\t')

def download_TCGA_phenotype_data(cancer_type_selector):
    """Saves TCGA phenotype data for the specifies @cancer_type obtained from UCSC Xena in /data.

    Parameters
    ----------
    cancer_type_selector : CancerTypeSelector
        Specifies the cohort that the phenotype file is to be downloaded for.
    """
    cwd = os.getcwd()
    url = ""
    if cancer_type_selector == CancerTypeSelector.BLCA:
        url = "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-BLCA.GDC_phenotype.tsv.gz"
    pheno_data = pd.read_csv(url, delimiter='\t')
    pheno_data =  pheno_data[pheno_data['sample_type.samples'] == 'Primary Tumor']
    pheno_data.to_csv(os.path.join(cwd, 'data', 'TCGA-'+str(cancer_type_selector)+'.GDC_phenotype.tsv'), sep='\t')

def download_known_tfs():
    """Saves known human transcription factors obtained from humantfs.ccbr.utoronto.ca in /data.
    """
    cwd = os.getcwd()
    df = pd.read_csv('http://humantfs.ccbr.utoronto.ca/download/v_1.01/TFs_Ensembl_v_1.01.txt', delimiter='\t', index_col=0)
    df.to_csv(os.path.join(cwd, 'data', 'regulators.csv'), sep='\t')

def get_expression_data(cancer_type_selector):
    """Loads the expression data for the selected cancer type.

    Parameters
    ----------
    cancer_type_selector : CancerTypeSelector
        Specifies for which cancer type the phenotypes should be loaded.
        
    Returns
    -------
    expression_data : pd.DataFrame
        Expression data (indices are sample IDs, column names are gene IDs).
    """
    cwd = os.getcwd()
    try:
        expression_data = pd.read_csv(os.path.join(cwd, 'data', 'TCGA-'+str(cancer_type_selector)+'.htseq_fpkm.tsv'), sep='\t', header=0, index_col=0)
    except FileNotFoundError:
        download_TCGA_expression_data(cancer_type_selector)
        expression_data = pd.read_csv(os.path.join(cwd, 'data', 'TCGA-'+str(cancer_type_selector)+'.htseq_fpkm.tsv'), sep='\t', header=0, index_col=0)

    return expression_data # TODO check if it really works if we do not filter on Primary Tumor here, but only in samples

def get_pheno_data(cancer_type_selector):
    """Loads the phenotype data for the selected cancer type.

    Parameters
    ----------
    cancer_type_selector : CancerTypeSelector
        Specifies for which cancer type the phenotypes should be loaded.

    Returns
    -------
    pheno_data : pd.DataFrame
        Expression data (indices are sample IDs, column names are gene IDs).
    """
    cwd = os.getcwd()
    try:
        pheno_data = pd.read_csv(os.path.join(cwd, 'data', 'TCGA-'+str(cancer_type_selector)+'.GDC_phenotype.tsv'), sep='\t', header=0, index_col=0,
        dtype = {'gender.demographic': str,'race.demographic': str, 'age_at_initial_pathologic_diagnosis': int, 'submitter_id.samples': str})
    except FileNotFoundError:
        download_TCGA_phenotype_data(cancer_type_selector)
        pheno_data = pd.read_csv(os.path.join(cwd, 'data', 'TCGA-'+str(cancer_type_selector)+'.GDC_phenotype.tsv'), sep='\t', header=0, index_col=0,
        dtype = {'gender.demographic': str,'race.demographic': str, 'age_at_initial_pathologic_diagnosis': int, 'submitter_id.samples': str})
    return pheno_data # TODO check if it really works if we do not filter on Primary Tumor here, but only in samples


def get_conf_partition(pheno_data, confounder_selector):
    """Returns two lists with the first containing string-identifiers for the blocks of the requested confounder 
        and the second containing the sample ids corresponding to the blocks.

        Parameters
        ----------
        pheno_data : pd.DataFrame
            Data frame containing phenotypic information. One row per sample, one column per attribute.

        confounder_selector : ConfounderSelector
            Confounder attribute that is to be used to build the partition.

        Returns
        -------
        conf_partition : list
            Contains the blocks belonging to the confounder-based partition.
    """
    indices = None
    if confounder_selector == ConfounderSelector.SEX:
        female_samples = pheno_data.loc[pheno_data['gender.demographic'] == 'female']['submitter_id.samples']
        male_samples = pheno_data.loc[pheno_data['gender.demographic'] == 'male']['submitter_id.samples']
        blocks, conf_partition = ['female', 'male'], [female_samples.tolist(), male_samples.tolist()]
    if confounder_selector == ConfounderSelector.RACE:
        asian_samples = pheno_data.loc[pheno_data['race.demographic'] == 'asian']['submitter_id.samples']
        african_samples = pheno_data.loc[pheno_data['race.demographic'] == 'black or african american']['submitter_id.samples']
        white_samples = pheno_data.loc[pheno_data['race.demographic'] == 'white']['submitter_id.samples']
        blocks, conf_partition = ['asian', 'african', 'white'], [asian_samples.tolist(), african_samples.tolist(), white_samples.tolist()]
    if confounder_selector == ConfounderSelector.AGE:
        low_age_samples = pheno_data.loc[pheno_data['age_at_initial_pathologic_diagnosis'] < 50]['submitter_id.samples']
        high_age_samples = pheno_data.loc[pheno_data['age_at_initial_pathologic_diagnosis'] > 70]['submitter_id.samples']
        blocks, conf_partition = ['low_age', 'high_age'], [low_age_samples.tolist(), high_age_samples.tolist()]
    return conf_partition

def get_n_random_partitions(n, samples, conf_partition):
    """Returns n random partitions each containing blocks of the same size as the corresponding blocks in the
    confounder based partition.

    Parameters
    ----------
    n : int
        Specifies the number of random partitions that should be generated.

    samples : pd.DataFrame
        Contains all sample identifiers.

    conf_partition : list
        List of blocks as pd.DataFrames with one column containing the sample identifiers belonging to the block.
        
    Returns
    -------
    partitions : list
        List of random partitions.
    """
    partitions = []
    samples_cpy = samples.copy()
    for k in range(n):
        cur = []
        for i in range(len(conf_partition)):
            block = samples_cpy.sample(n=len(conf_partition[i]), replace=True)
            cur.append(block)
        partitions.append(cur)
    return partitions