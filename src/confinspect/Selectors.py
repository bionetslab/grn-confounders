from enum import Enum
from .GENIE3Wrapper import GENIE3Wrapper
from .ARACNEWrapper import ARACNEWrapper
from .WGCNAWrapper import WGCNAWrapper
from .CEMiWrapper import CEMiWrapper
from .GRNBOOST2Wrapper import GRNBOOST2Wrapper
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
    ALL = 'ALL'

    def __str__(self):
        return self.value

class AlgorithmSelector(Enum):
    """Enum specifying which network inference algorithm should be used."""
    ARACNE = 'ARACNE'
    GENIE3 = 'GENIE3'
    WGCNA = 'WGCNA'
    CEMITOOL = 'CEMITOOL'
    GRNBOOST2 = 'GRNBOOST2'
    
    def __str__(self):
        return self.value

class BlockType(Enum):
    """Enum listing possible block types of confounders."""
    QUARTILE = 'QUARTILE'
    CATEGORY = 'CATEGORY'
    ALL = 'ALL'
    
    def __str__(self):
        return self.value

class Role(Enum):
    """Enum listing possible role identifiers of fields, i.e. confounders or variables."""
    CONFOUNDER = 'CONFOUNDER'
    VARIABLE = 'VARIABLE'
    
    def __str__(self):
        return self.value

def get_algorithm_wrapper(algorithm_selector):
    """Returns the appropriate algorithm based on the selection.

    Parameters
    ----------
    algorithm_selector : AlgorithmSelector
        Specifies which algorithm should be used.
    cwd : str
        Current working directory containing all input and output directories.
    """
    if algorithm_selector == AlgorithmSelector.GENIE3:
        return GENIE3Wrapper()
    elif algorithm_selector == AlgorithmSelector.ARACNE:
        return ARACNEWrapper()
    elif algorithm_selector == AlgorithmSelector.WGCNA:
        return WGCNAWrapper()
    elif algorithm_selector == AlgorithmSelector.CEMITOOL:
        return CEMiWrapper()
    elif algorithm_selector == AlgorithmSelector.GRNBOOST2:
        return GRNBOOST2Wrapper()
