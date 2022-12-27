"""
def _new_str(self):
    return self.value

def dynamicConfounderSelector(values):
    
    dynamic_enum = Enum('ConfounderSelector', values)
    dynamic_enum.__str__ = funcType(_new_str, dynamic_enum, dynamicConfounderSelector)
    return dynamic_enum

def dynamicCancerTypeSelector(values):
    #Dynamic enum specifying the cohorts to be examined. Is initialized based on the user input.
    #Parameters
    #----------
    #values : dict
        #Dictionary of enum elements with string values.

    #Returns
    #----------
    #dynamic_enum : Enum
        #Custom enum with elements that represent the cohorts to be examined.
    
    dynamic_enum = Enum('CancerTypeSelector', values)
    dynamic_enum.__str__ = funcType(_new_str, dynamic_enum, dynamicCancerTypeSelector)
    return dynamic_enum

def parse_input_file():
    args = {}
    args['combine'] = False
    f = open(os.path.join(os.getcwd(), "input_parameters.txt"))
    input_data = f.readlines()
    for line in input_data:
        key, value = line.split(":")
        if key == 'combine'.strip():
            args[key.strip()] = True
        elif key.strip() in ['N_from', 'N_to', 'M_from', 'M_to', 'k']:
            args[key.strip()] = int(value.strip())
        elif key == 'mode'.strip():
            assert value.strip() == 'par' or value.strip() == 'seq'
            args[key.strip()] = value.strip()
        else:
            args[key.strip()] = value.strip()
    f.close()
    assert all([el in args.keys() for el in ['ct', 'conf', 'alg', 'N_from', 'N_to', 'M_from', 'M_to', 'k', 'mode']])
    args = access_with_dot(args)
    return args

class access_with_dot(dict):
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

"""

def download_TCGA_expression_data(cancer_type_selector):
    """Saves TCGA gene expression RNAseq - HTSeq - FPKM data for the specifies @cancer_type obtained from UCSC Xena in /data.

    Parameters
    ----------
    cancer_type_selector : str
        Specifies the cohort that the phenotype file is to be downloaded for.
    """
    cwd = os.getcwd()
    url = ""
    if cancer_type_selector in [str(val) for val in list(TCGACancerTypeSelector)]:
        url = f'https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-{str(cancer_type_selector)}.htseq_fpkm.tsv.gz'
        df = pd.read_csv(url, delimiter='\t', index_col='Ensembl_ID').T
        df.to_csv(os.path.join(cwd, 'data', 'TCGA-'+str(cancer_type_selector)+'.htseq_fpkm.tsv'), sep='\t')

def download_TCGA_phenotype_data(cancer_type_selector):
    """Saves TCGA phenotype data for the specifies @cancer_type obtained from UCSC Xena in /data.

    Parameters
    ----------
    cancer_type_selector : str
        Specifies the cohort that the phenotype file is to be downloaded for.
    """
    cwd = os.getcwd()
    url = ""
    if cancer_type_selector in [str(val) for val in list(TCGACancerTypeSelector)]:
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
