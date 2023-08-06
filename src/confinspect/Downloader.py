import pandas as pd
import os
import argparse
from . import Selectors

def download_TCGA_expression_data(cancer_type_selector):
    """Saves TCGA gene expression RNAseq - HTSeq - FPKM data for the specifies @cancer_type obtained from UCSC Xena in /data.

    Parameters
    ----------
    cancer_type_selector : str
        Specifies the cohort that the phenotype file is to be downloaded for.
    """
    cwd = os.getcwd()
    if not os.path.exists(os.path.join(cwd, 'data')):
        os.mkdir(os.path.join(cwd, 'data'))
    url = ""
    if cancer_type_selector in [str(val) for val in list(Selectors.TCGACancerTypeSelector)]:
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
    if cancer_type_selector in [str(val) for val in list(Selectors.TCGACancerTypeSelector)]:
        url = f'https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-{str(cancer_type_selector)}.GDC_phenotype.tsv.gz'
        pheno_data = pd.read_csv(url, delimiter='\t', index_col='submitter_id.samples')
        pheno_data =  pheno_data[pheno_data['sample_type.samples'] == 'Primary Tumor']
        
        if not os.path.exists(os.path.join(cwd, 'data')):
            os.mkdir(os.path.join(cwd, 'data'))
        pheno_data.to_csv(os.path.join(cwd, 'data', 'TCGA-'+str(cancer_type_selector)+'.GDC_phenotype.tsv'), sep='\t')

def download_known_tfs():
    """Saves known human transcription factors obtained from humantfs.ccbr.utoronto.ca in /data.
    """
    cwd = os.getcwd()
    if not os.path.exists(os.path.join(cwd, 'data')):
        os.mkdir(os.path.join(cwd, 'data'))
    df = pd.read_csv('http://humantfs.ccbr.utoronto.ca/download/v_1.01/TFs_Ensembl_v_1.01.txt', delimiter='\t', index_col=0)
    df.to_csv(os.path.join(cwd, 'data', 'regulators_.csv'), sep=',')

def download_protein_coding_genes():
    """Saves protein coding genes obtained from http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/locus_groups/protein-coding_gene.txt
    in /data."""
    cwd = os.getcwd()
    if not os.path.exists(os.path.join(cwd, 'data')):
        os.mkdir(os.path.join(cwd, 'data'))
    df = pd.read_csv('http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/locus_groups/protein-coding_gene.txt', delimiter='\t', index_col=False, dtype=str)
    df = df[['ensembl_gene_id', 'symbol']]
    df.to_csv(os.path.join(cwd, 'data', 'protein-coding_gene.csv'), sep=',', index=False)

