from enum import Enum
import GENIE3Wrapper
import ARACNEWrapper
import LSCONWrapper
import pandas as pd
import os

class AlgorithmSelector(Enum):
    """Enum specifying which network inference algorithm should be used."""
    GENIE3 = 'GENIE3'
    ARACNe = 'ARACNe'
    LSCON = 'LSCON'

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
    elif algorithm_selector == AlgorithmSelector.ARACNe:
        return ARACNe()
    elif algorithm_selector == AlgorithmSelector.LSCON:
        return LSCON()

def download_TCGA_expression_data(cancer_type_selector):
    """Saves TCGA gene expression RNAseq - HTSeq - FPKM data for the specifies @cancer_type obtained from UCSC Xena in /data.
    Parameters
    ----------
    cancer_type_selector : CancerTypeSelector
        Specifies the cohort that the phenotype file is to be downloaded for.
    """
    cwd = os.path.join(os.path.dirname(__file__))
    url = ""
    if cancer_type_selector == CancerTypeSelector.BLCA:
        url = "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-BLCA.htseq_fpkm.tsv.gz"
    df = pd.read_csv(url, delimiter='\t', index_col='Ensembl_ID').T
    df.to_csv(os.path.join(cwd, '..', 'data', 'TCGA-'+str(cancer_type_selector)+'.htseq_fpkm.tsv'), sep='\t')

def download_TCGA_phenotype_data(cancer_type_selector):
    """Saves TCGA phenotype data for the specifies @cancer_type obtained from UCSC Xena in /data.
    Parameters
    ----------
    cancer_type_selector : CancerTypeSelector
        Specifies the cohort that the phenotype file is to be downloaded for.
    """
    cwd = os.path.join(os.path.dirname(__file__))
    url = ""
    if cancer_type_selector == CancerTypeSelector.BLCA:
        url = "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-BLCA.GDC_phenotype.tsv.gz"
    pheno_data = pd.read_csv(url, delimiter='\t')
    pheno_data =  pheno_data[pheno_data['sample_type.samples'] == 'Primary Tumor']
    pheno_data.to_csv(os.path.join(cwd, '..', 'data', 'TCGA-'+str(cancer_type_selector)+'.GDC_phenotype.tsv'), sep='\t')

def download_known_tfs():
    """Saves known human transcription factors obtained from humantfs.ccbr.utoronto.ca in /data.
    """
    df = pd.read_csv('http://humantfs.ccbr.utoronto.ca/download/v_1.01/TFs_Ensembl_v_1.01.txt', delimiter='\t', index_col=0)
    df.to_csv(os.path.join(cwd, '..', 'data', 'regulators.csv'), sep='\t')

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
    expression_data = pd.read_csv(os.path.join(cwd, '..', 'data', 'TCGA-'+str(cancer_type_selector)+'.htseq_fpkm.tsv'), sep='\t')
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
    pheno_data = pd.read_csv(os.path.join(cwd, '..', 'data', 'TCGA-'+str(cancer_type_selector)+'.GDC_phenotype.tsv'), sep='\t')
    return pheno_data # TODO check if it really works if we do not filter on Primary Tumor here, but only in samples


def get_conf_partition(pheno_data, confounder):
    """Returns two np.arrays with the first containing string-identifiers for the blocks of the requested confounder 
        and the second containing the sample ids corresponding to the blocks.
            
        Parameters
        ----------
        pheno_data : pd.DataFrame
            Data frame containing phenotypic information. One row per sample, one column per attribute.

        confounder : str
            Confounder attribute that is to be used to build the partition.
        
        Returns
        -------
        blocks : np.array
            At index i, contains the str-identifier of the block at index i in the @conf_partition.
        
        conf_partition : np.array
            Contains the blocks belonging to the confounder-based partition.
    """
    indices = None

    # split set of sample ids based on confounder expression
    if confounder == 'sex':
        gender_col_index = np.where(pheno_data[0, :] == 'gender.demographic')[0][0]
        female_data = pheno_data[pheno_data[:, gender_col_index] == 'female']
        male_data = pheno_data[pheno_data[:, gender_col_index] == 'male']

        female_samples = female_data[:, 0]
        male_samples = male_data[:, 0]
        blocks, conf_partition = ['female', 'male'], [female_samples, male_samples]

    if confounder == 'ethnicity':
        ethn_col_index = np.where(pheno_data[0, :] == 'ethnicity.demographic')[0][0]
        hisp_lat_data = pheno_data[pheno_data[:, ethn_col_index] == 'hispanic or latino']
        nonHisp_nonLat_data = pheno_data[pheno_data[:, ethn_col_index] == 'not hispanic or latino']

        hisp_lat_samples = hisp_lat_data[:, 0]
        nonHisp_nonLat_samples = nonHisp_nonLat_data[:, 0]
        blocks, conf_partition = ['hisp_lat', 'non_hisp_lat'], [hisp_lat_samples, nonHisp_nonLat_samples]

    if confounder == 'race':
        race_col_index = np.where(pheno_data[0, :] == 'race.demographic')[0][0]
        asian_data = pheno_data[pheno_data[:, race_col_index] == 'asian']
        african_data = pheno_data[pheno_data[:, race_col_index] == 'black or african american']
        white_data = pheno_data[pheno_data[:, race_col_index] == 'white']

        asian_samples = asian_data[:, 0]
        african_samples = african_data[:, 0]
        white_samples = white_data[:, 0]
        blocks, conf_partition = ['asian', 'african', 'white'], [asian_samples, african_samples, white_samples]

    return conf_partition

def get_n_random_partitions(n, samples, conf_partition):
    partitions = []
    pos = 0
    cuts = []
    for block in self.conf_partition:
        cuts.append(pos+len(block))
        pos = pos + len(block) # this might seem weird, but the np.split function is very handy and it needs the indices where the cuts should be placed
    for k in range(n): # TODO
        np.random.shuffle(samples.copy())
        partitions.append(np.split(samples_cpy, cuts[:-1]))
    return partitions