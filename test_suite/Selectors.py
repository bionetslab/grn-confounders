from enum import Enum
"""
from .GENIE3Wrapper import GENIE3Wrapper
from .ARACNEWrapper import ARACNEWrapper
"""
from WGCNAWrapper import WGCNAWrapper
"""
from .CEMiWrapper import CEMiWrapper
from .GRNBOOST2Wrapper import GRNBOOST2Wrapper
from .sdcorGCNWrapper import sdcorGCNWrapper
"""
import pandas as pd
import numpy as np
import os

class TCGACancerTypeSelector(Enum):
    """Enum listing TCGA cancer type selectors by study abbreviations."""

    ACC = 'ACC'
    LAML = 'LAML'
    CHOL = 'CHOL'
    BLCA = 'BLCA'
    UCEC = 'UCEC'
    ESCA = 'ESCA'
    KICH = 'KICH'
    DLBC = 'DLBC'
    LIHC = 'LIHC'
    LGG = 'LGG'
    LUAD = 'LUAD'
    SKCM = 'SKCM'
    MESO = 'MESO'
    UVM = 'UVM'
    OV = 'OV'
    PAAD = 'PAAD'
    PRAD = 'PRAD'
    SARC = 'SARC'
    TGCT = 'TGCT'
    THYM = 'THYM'
    THCA = 'THCA'
    UCS = 'UCS'
    PCPG = 'PCPG'
    GBM = 'GBM'
    COAD = 'COAD'
    STAD = 'STAD'
    READ = 'READ'
    BRCA = 'BRCA'
    LUSC = 'LUSC'
    HNSC = 'HNSC'
    CESC = 'CESC'
    KIRP = 'KIRP'
    KIRC = 'KIRC'

    def __str__(self):
        return self.value

class ConfounderSelector(Enum):
    """Enum listing predefined confounder selectors to avoid confusing custom string confounders and predefined confounders."""

    TYPE = 'TYPE'
    NONE = 'NONE'

    def __str__(self):
        return self.value

class AlgorithmSelector(Enum):
    """Enum specifying which network inference algorithm should be used."""
    ARACNE = 'ARACNE'
    GENIE3 = 'GENIE3'
    WGCNA = 'WGCNA'
    CEMI = 'CEMI'
    GRNBOOST2 = 'GRNBOOST2'
    SDCORGCN = 'SDCORGCN'
    CUSTOMGRN = 'CUSTOMGRN'
    CUSTOMGCN = 'CUSTOMGCN'
    
    def __str__(self):
        return self.value

class BlockType(Enum):
    """Enum listing possible block types of confounders."""
    QUARTILE = 'QUARTILE'
    CATEGORY = 'CATEGORY'
    All = 'ALL'
    
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
    elif algorithm_selector == AlgorithmSelector.SDCORGCN:
        return sdcorGCNWrapper()
    elif algorithm_selector == AlgorithmSelector.CUSTOMGRN:
        return CUSTOMGRNWrapper()
    elif algorithm_selector == AlgorithmSelector.CUSTOMGCN:
        return CUSTOMGCNWrapper()

def get_expression_data(cancer_type_selector, ged, sep='\t'):
    """Loads the expression data for the selected cancer type. Only leaves columns of protein-codeing genes, provided
    in the file protein-coding_gene.csv, in expression_data. Removes gene version identifiers from gene ensembl IDs.

    Parameters
    ----------
    cancer_type_selector : str
        Specifies for which cancer type the phenotypes should be loaded.

    ged : str
        String describing file name of gene expression data file in data directory. 
        
    Returns
    -------
    expression_data : pd.DataFrame
        Expression data (indices are sample IDs, column names are gene IDs).
    """
    cwd = os.getcwd()
    expression_data = pd.read_csv(os.path.join(cwd, 'data', ged), sep=sep, header=0, index_col=0)
    print('Remove version identifiers from gene symbols in expression data for cohort ' + str(cancer_type_selector) + '...')
    expression_data.columns = expression_data.columns.str.split('.').str[0].tolist()
    print('Only leave protein-coding genes in expression data set for cohort ' + str(cancer_type_selector) + '...')
    pcg = pd.read_csv(os.path.join('data', 'protein-coding_gene.csv'))
    mask = expression_data.columns.intersection(pcg['ensembl_gene_id'].values)
    expression_data = expression_data[mask]
    return expression_data

def get_pheno_data(cancer_type_selector, pt, sep='\t'):
    """Loads the phenotype data for the selected cancer type and adds a column containing cancer_type_selector. Only leave
    rows (samples) in pheno_data that have non-NaN values in the columns gender.demographic, race.demographic, 
    age_at_initial_pathologic_diagnosis and tumor_stage.diagnoses.

    Parameters
    ----------
    cancer_type_selector : str
        Specifies for which cancer type the phenotypes should be loaded.

    Returns
    -------
    pheno_data : pd.DataFrame
        Expression data (indices are sample IDs, column names are gene IDs).
    """
    cwd = os.getcwd()
    pheno_data = pd.read_csv(os.path.join(cwd, 'data', pt), sep=sep, header=0, index_col=0)
    assert len(pheno_data.iloc[0]) == len(pheno_data.iloc[0].values)
    pheno_data = pheno_data.set_index(pheno_data.columns[0])
    pheno_data['cohort'] = str(cancer_type_selector)
    print('Filter Primary Tumor samples in pheno data for cohort ' + str(cancer_type_selector) + '...')
    pheno_data =  pheno_data[pheno_data['sample_type.samples'] == 'Primary Tumor']
    return pheno_data

def get_conf_partition(pheno_data_orig, block_type, pheno_field, rank=0, min_block_size=20):
    """Returns two lists with the first containing string-identifiers for the blocks of the requested confounder 
    and the second containing the sample ids corresponding to the blocks. For the age confounder, the lower and upper quartiles are
    computed separately per cohort and are then combined into the resulting upper and lower fragments.

    Parameters
    ----------
    pheno_data_orig : pd.DataFrame
        Data frame containing phenotypic information. One row per sample, one column per attribute.

    pheno_field : str
        Field in pheno type file to be used to induce the partition.

    block_type : BlockType
        QUARTILE or CATEGORY; defines how to create the partition.

    Returns
    -------
    conf_partition : list
        Contains the blocks belonging to the confounder-based partition.
    """
    pheno_data = pheno_data_orig.copy()
    indices = None
    blocks = []
    conf_partition = []
    pheno_data = pheno_data[pheno_data[pheno_field] != 'not reported']
    pheno_data = pheno_data[pheno_data[pheno_field].notna()]
    if pheno_field == 'tumor_stage.diagnoses':
        pheno_data = pheno_data[pheno_data[pheno_field] != 'stage x']
        pheno_data.loc[pheno_data['tumor_stage.diagnoses'].str.strip().isin(['stage ia', 'stage ib', 'stage ic']), pheno_field] = 'stage i'
        pheno_data.loc[pheno_data['tumor_stage.diagnoses'].str.strip().isin(['stage iia', 'stage iib', 'stage iic']), pheno_field] = 'stage ii'
        pheno_data.loc[pheno_data['tumor_stage.diagnoses'].str.strip().isin(['stage iiia', 'stage iiib', 'stage iiic', 'stage iv', 'stage iva', 'stage ivb', 'stage ivc']), pheno_field] = 'stage iii'
    if pheno_field == ConfounderSelector.NONE:
        samples = pheno_data.index.tolist()
        conf_partition.append(('all', samples))
    elif block_type == BlockType.CATEGORY:
        blocks = sorted(list(set(pheno_data[pheno_field].str.strip().values)))
        for block_attr in blocks:
            samples = pheno_data.loc[pheno_data[pheno_field].str.strip() == block_attr].index.tolist()
            if len(samples) >= min_block_size:
                conf_partition.append((block_attr, samples))
    elif block_type == BlockType.QUARTILE:
        samples_lower = []
        samples_upper = []
        for cohort in set(pheno_data['cohort'].str.strip().values):
            pheno_cohort = pheno_data[pheno_data['cohort'] == cohort]
            lower, upper = pheno_cohort[pheno_field].quantile(0.25), pheno_cohort[pheno_field].quantile(0.75)
            samples_lower.extend(pheno_cohort.loc[pheno_cohort[pheno_field] <= lower].index.tolist())
            samples_upper.extend(pheno_cohort.loc[pheno_cohort[pheno_field] > upper].index.tolist())
        if len(samples_lower) >= min_block_size and len(samples_upper) >= min_block_size:
            conf_partition.append(('lower', samples_lower))
            conf_partition.append(('upper', samples_upper))
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
    if conf_sel == ConfounderSelector.NONE:
        partitions.append(samples)
        return partitions
    for k in range(n_from, n_to):
        samples_cpy = pd.DataFrame(samples.copy())
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
    return partitions

def align_data(expression_datasets, pheno_datasets):
    """Data alignment. Remove such samples from the expression_data files that are not in the pheno_data files and vice versa.
    Remove all genes where standard deviation is 0."""
    for sel in expression_datasets.keys():
        print('Align expression data and phenotype data on samples for cohort ' + str(sel) + '...')
        pheno_datasets[sel] = pheno_datasets[sel][pheno_datasets[sel].index.isin(expression_datasets[sel].index)]
        samples = pheno_datasets[sel].index.values        
        expression_datasets[sel] = expression_datasets[sel].loc[samples]
        print('Remove genes where standard deviation of expression data is 0 for cohort ' + str(sel) + '...')
        expression_datasets[sel] = expression_datasets[sel].loc[:, (expression_datasets[sel].std() != 0)]
        pheno_datasets[sel] = pheno_datasets[sel][pheno_datasets[sel].index.isin(expression_datasets[sel].index)]
    return expression_datasets, pheno_datasets