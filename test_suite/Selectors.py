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

def get_algorithm_wrapper(algorithm_selector):
    """Returns the appropriate algorithm based on the selection.
    Parameters
    ----------
    algorithm_selector : AlgorithmSelector
        Specifies which algorithm should be used.
    """
    if algorithm_selector == AlgorithmSelector.GENIE3:
        return GENIEWrapper()
    elif algorithm_selector == AlgorithmSelector.ARACNe:
        return ARACNe()
    elif algorithm_selector == AlgorithmSelector.LSCON:
        return LSCON()

def download_TCGA_expression_data(cancer_type_selector):
    """Returns TCGA gene expression RNAseq - HTSeq - FPKM data for the specifies @cancer_type obtained from USCS Xena.
    Parameters
    ----------
    cancer_type_selector : CancerTypeSelector
        Specifies which algorithm should be used.

    Returns
    -------
    raw_expression_data : pd.DataFrame
        A data frame with gene symbols as indices and column names whose entries correspond to
        non-normalized gene expression data.
    """
    cwd = os.path.join(os.path.dirname(__file__))
    url = ""
    if cancer_type_selector == CancerTypeSelector.BLCA:
        url = "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-BLCA.htseq_fpkm.tsv.gz"
    df = pd.read_csv(url, delimiter='\t', index_col='Ensembl_ID').T
    df.to_csv(os.path.join(cwd, '..', 'data', 'TCGA-BLCA.htseq_fpkm.tsv'), sep='\t')
