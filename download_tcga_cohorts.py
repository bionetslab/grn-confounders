import pandas as pd
import numpy as np
import os
from test_suite import Selectors
import argparse

def download_TCGA_expression_data(cancer_type_selector):
    """Saves TCGA gene expression RNAseq - HTSeq - FPKM data for the specifies @cancer_type obtained from UCSC Xena in /data.

    Parameters
    ----------
    cancer_type_selector : str
        Specifies the cohort that the phenotype file is to be downloaded for.
    """
    cwd = os.getcwd()
    url = ""
    if cancer_type_selector in [str(val) for val in list(Selectors.TCGACancerTypeSelector)]:
        url = f'https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-{str(cancer_type_selector)}.htseq_fpkm.tsv.gz'
        df = pd.read_csv(url, delimiter='\t', index_col='Ensembl_ID').T
        df.to_csv(os.path.join(cwd, 'data', 'TCGA-'+str(cancer_type_selector)+'.htseq_fpkm.tsv'), sep=',')

def download_TCGA_phenotype_data(cancer_type_selector):
    """Saves TCGA phenotype data for the specifies @cancer_type obtained from UCSC Xena in /data.

    Parameters
    ----------
    cancer_type_selector : str
        Specifies the cohort that the phenotype file is to be downloaded for.
    """
    cwd = os.getcwd()
    url = ""
    if cancer_type_selector in [str(val) for val in list(Selectors.TCGACancerTypeSelector)]:
        url = f'https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-{str(cancer_type_selector)}.GDC_phenotype.tsv.gz'
        pheno_data = pd.read_csv(url, delimiter='\t', index_col='submitter_id.samples')
        pheno_data =  pheno_data[pheno_data['sample_type.samples'] == 'Primary Tumor']
        pheno_data.to_csv(os.path.join(cwd, 'data', 'TCGA-'+str(cancer_type_selector)+'.GDC_phenotype.tsv'), sep=',')

def download_known_tfs():
    """Saves known human transcription factors obtained from humantfs.ccbr.utoronto.ca in /data.
    """
    cwd = os.getcwd()
    df = pd.read_csv('http://humantfs.ccbr.utoronto.ca/download/v_1.01/TFs_Ensembl_v_1.01.txt', delimiter='\t', index_col=0)
    df.to_csv(os.path.join(cwd, 'data', 'regulators.csv'), sep=',')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
                    prog = 'download_TCGA_cohorts',
                    description = 'Download TCGA gene expression data and phenotype data and list of knon transcription factors as list of regulators for GRN methods.')    
    parser.add_argument('-tfs', action='store_true', help='If set, download list of human known transcription factors')
    parser.add_argument('-ct', required=True, nargs='+', help='Use TCGA study abbreviations to specify the cohorts to be downloaded')
    args = parser.parse_args()
    for ct_sel in list(args.ct):
        download_TCGA_expression_data(ct_sel)
        download_TCGA_phenotype_data(ct_sel)
    if args.tfs:
        download_known_tfs()
