from enum import Enum
from .GENIE3Wrapper import GENIE3Wrapper
from .ARACNEWrapper import ARACNEWrapper
from .WGCNAWrapper import WGCNAWrapper
from .CEMiWrapper import CEMiWrapper
from .GRNBOOST2Wrapper import GRNBOOST2Wrapper
import pandas as pd
import numpy as np
import os

class AlgorithmSelector(Enum):
    """Enum specifying which network inference algorithm should be used."""
    ARACNE = 'ARACNE'
    GENIE3 = 'GENIE3'
    WGCNA = 'WGCNA'
    CEMI = 'CEMI'
    GRNBOOST2 = 'GRNBOOST2'
    
    def __str__(self):
        return self.value

class CancerTypeSelector(Enum):
    """Enum specifying which cancer type should be investigated."""
    #BLCA = 'BLCA'
    PCPG = 'PCPG'
    GBM = 'GBM'
    COAD = 'COAD'
    STAD = 'STAD'
    READ = 'READ'
    BRCA = 'BRCA'
    #LUAD = 'LUAD'
    LUSC = 'LUSC'
    HNSC = 'HNSC'
    CESC = 'CESC'
    KIRP = 'KIRP'
    KIRC = 'KIRC'
    CHOL = 'CHOL'
    LIHC = 'LIHC'
    KICH = 'KICH'
    ACC = 'ACC'
    #PRAD = 'PRAD'
    #SKCM = 'SKCM'

    def __str__(self):
        return self.value

class ConfounderSelector(Enum):
    """Enum specifying the confounder whose effect is to be examined."""
    SEX = 'sex'
    RACE = 'race'
    AGE = 'age'
    STAGE = 'stage'
    TYPE = 'type'

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
    elif algorithm_selector == AlgorithmSelector.WGCNA:
        return WGCNAWrapper()
    elif algorithm_selector == AlgorithmSelector.CEMI:
        return CEMiWrapper()
    elif algorithm_selector == AlgorithmSelector.GRNBOOST2:
        return GRNBOOST2Wrapper()

def download_TCGA_expression_data(cancer_type_selector):
    """Saves TCGA gene expression RNAseq - HTSeq - FPKM data for the specifies @cancer_type obtained from UCSC Xena in /data.

    Parameters
    ----------
    cancer_type_selector : CancerTypeSelector
        Specifies the cohort that the phenotype file is to be downloaded for.
    """
    cwd = os.getcwd()
    url = ""
    if cancer_type_selector in list(CancerTypeSelector):
        url = f'https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-{str(cancer_type_selector)}.htseq_fpkm.tsv.gz'
        df = pd.read_csv(url, delimiter='\t', index_col='Ensembl_ID').T
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
    if cancer_type_selector in list(CancerTypeSelector):
        url = f'https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-{str(cancer_type_selector)}.GDC_phenotype.tsv.gz'
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
    """Loads the expression data for the selected cancer type. Only leaves columns of protein-codeing genes, provided
    in the file protein-coding_gene.csv, in expression_data. Removes gene version identifiers from gene ensembl IDs.

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
    print('Remove version identifiers from gene symbols in expression data for cohort ' + str(cancer_type_selector) + '...')
    expression_data.columns = expression_data.columns.str.split('.').str[0].tolist()
    print('Only leave protein-coding genes in expression data set for cohort ' + str(cancer_type_selector) + '...')
    pcg = pd.read_csv(os.path.join('data', 'protein-coding_gene.csv'))
    mask = expression_data.columns.intersection(pcg['ensembl_gene_id'].values)
    expression_data = expression_data[mask]
    return expression_data

def get_pheno_data(cancer_type_selector):
    """Loads the phenotype data for the selected cancer type and adds a column containing cancer_type_selector. Only leave
    rows (samples) in pheno_data that have non-NaN values in the columns gender.demographic, race.demographic, 
    age_at_initial_pathologic_diagnosis and tumor_stage.diagnoses.

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
        dtype = {'gender.demographic': str,'race.demographic': str, 'age_at_initial_pathologic_diagnosis': float, 'submitter_id.samples': str})
    except FileNotFoundError:
        download_TCGA_phenotype_data(cancer_type_selector)
        pheno_data = pd.read_csv(os.path.join(cwd, 'data', 'TCGA-'+str(cancer_type_selector)+'.GDC_phenotype.tsv'), sep='\t', header=0, index_col=0,
        dtype = {'gender.demographic': str,'race.demographic': str, 'age_at_initial_pathologic_diagnosis': float, 'submitter_id.samples': str})
    pheno_data['cohort'] = str(cancer_type_selector)
    print('Filter Primary Tumor samples in pheno data for cohort ' + str(cancer_type_selector) + '...')
    pheno_data =  pheno_data[pheno_data['sample_type.samples'] == 'Primary Tumor']

    return pheno_data

def get_conf_partition(pheno_data_orig, confounder_selector):
    """Returns two lists with the first containing string-identifiers for the blocks of the requested confounder 
    and the second containing the sample ids corresponding to the blocks.

    Parameters
    ----------
    pheno_data_orig : pd.DataFrame
        Data frame containing phenotypic information. One row per sample, one column per attribute.

    confounder_selector : ConfounderSelector
        Confounder attribute that is to be used to build the partition.

    Returns
    -------
    conf_partition : list
        Contains the blocks belonging to the confounder-based partition.
    """
    pheno_data = pheno_data_orig.copy()
    indices = None
    blocks = []
    conf_partition = []
    pheno_field = ''  
    if confounder_selector == ConfounderSelector.SEX:
        pheno_field = 'gender.demographic'
    elif confounder_selector == ConfounderSelector.RACE:
        pheno_field = 'race.demographic'
    elif confounder_selector == ConfounderSelector.AGE:
        pheno_field = 'age_at_initial_pathologic_diagnosis'
    elif confounder_selector == ConfounderSelector.STAGE:
        pheno_field = 'tumor_stage.diagnoses'
        pheno_data = pheno_data[pheno_data[pheno_field] != 'stage x']
        pheno_data.loc[pheno_data['tumor_stage.diagnoses'].str.strip().isin(['stage ia', 'stage ib', 'stage ic'])] = 'stage i'
        pheno_data.loc[pheno_data['tumor_stage.diagnoses'].str.strip().isin(['stage iia', 'stage iib', 'stage iic'])] = 'stage ii'
        pheno_data.loc[pheno_data['tumor_stage.diagnoses'].str.strip().isin(['stage iiia', 'stage iiib', 'stage iiic', 'stage iv', 'stage iva', 'stage ivb', 'stage ivc'])] = 'stage iii'
    elif confounder_selector == ConfounderSelector.TYPE:
        pheno_field = 'cohort'
    pheno_data = pheno_data[pheno_data[pheno_field] != 'not reported']
    pheno_data = pheno_data[pheno_data[pheno_field].notna()]

    if confounder_selector != ConfounderSelector.AGE:
        blocks = list(set(pheno_data[pheno_field].str.strip().values))
        for block_attr in blocks:
            samples = pheno_data.loc[pheno_data[pheno_field].str.strip() == block_attr]['submitter_id.samples'].tolist()
            if len(samples) >= 20:
                conf_partition.append(samples)
            else:
                blocks.remove(block_attr)

    elif confounder_selector == ConfounderSelector.AGE:
        lower, upper = pheno_data[pheno_field].quantile(0.25), pheno_data[pheno_field].quantile(0.75)
        blocks = ['age_less_equal_'+str(lower), 'high_age_greater_'+str(upper)]
        samples_lower = pheno_data.loc[pheno_data[pheno_field] <= lower]['submitter_id.samples'].tolist()
        samples_upper = pheno_data.loc[pheno_data[pheno_field] > upper]['submitter_id.samples'].tolist()
        if len(samples_lower) >= 20 and len(samples_upper) >= 20:
            conf_partition.append(samples_lower)
            conf_partition.append(samples_upper)
        else:
            blocks = []

    with open('blocks_'+str(confounder_selector), 'a') as f:
        for i in range(len(blocks)):
            try:
                f.write(str(blocks[i])+': '+str(len(conf_partition[i]))+'\n')
            except IndexError:
                continue
    return conf_partition

def get_n_random_partitions(n_from, n_to, samples, conf_partition, ct_sel, conf_sel):
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
    partitions=[]
    for k in range(n_from, n_to):
        samples_cpy = samples.copy()
        cur = []
        try:
            part = pd.read_csv(os.path.join('partitions', f'rnd_part{k}_{ct_sel}_{conf_sel}'), header=None, index_col=False, dtype=str).values.tolist()
            begin = 0
            end = 0
            for i in range(len(conf_partition)):
                end += len(conf_partition[i])
                cur.append([item for sublist in part[begin:end] for item in sublist])
                begin += len(conf_partition[i])
        except FileNotFoundError:
            print(f'rnd_partition {k} not found. Create new partitions.')
            for i in range(len(conf_partition)):
                block = samples_cpy.sample(n=len(conf_partition[i]), replace=False)
                samples_cpy = samples_cpy[~samples_cpy.isin(block)]
                cur.append(block.values)
                block.to_csv(os.path.join('partitions', f'rnd_part{k}_{ct_sel}_{conf_sel}'), mode='a', header=False, index=False)
        partitions.append(cur)
        with open('block_specs', 'a') as f:
            f.write(str(ct_sel)+'_'+str(conf_sel)+'_'+str(len(samples_cpy))+'\n')
            for i in range(len(conf_partition)):
                f.write(str(len(cur[i]))+'\n')
            f.write('\n')
    return partitions
