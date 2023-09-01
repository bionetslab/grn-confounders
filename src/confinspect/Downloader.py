import pandas as pd
import numpy as np
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

def download_metabric_data():
    """Saves metabric expression data and phenotype data for metabric brca in /data.
    """
    cwd = os.getcwd()
    if not os.path.exists(os.path.join(cwd, 'data')):
        os.mkdir(os.path.join(cwd, 'data'))
    url = 'https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/brca_metabric/data_mrna_illumina_microarray.txt'
    df = pd.read_csv(url, delimiter='\t', index_col='Hugo_Symbol')

    # filter and translate gene identifiers to ensembl id
    df = df.drop('Entrez_Gene_Id', axis=1)
    df.index.name = None
    df = df.transpose()
    pcgs = pd.read_csv('http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/locus_groups/protein-coding_gene.txt', delimiter='\t', index_col=False, dtype=str)
    pcgs = pcgs[['ensembl_gene_id', 'symbol']]
    df = df[[col for col in pcgs['symbol'].values if col in df.columns]]
    df = df.loc[:,~df.columns.duplicated()].copy()
    assert pd.Series(df.columns).is_unique
    df = df.rename(columns=dict(zip(pcgs.symbol, pcgs.ensembl_gene_id)))
    df = df.dropna()
    # get pheno data of patients
    url = 'https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/brca_metabric/data_clinical_patient.txt'
    pheno = pd.read_csv(url, sep='\t', index_col=0)
    pheno = pheno.drop(['#Identifier to uniquely specify a patient.', '#STRING', '#1', 'PATIENT_ID'])
    # get pheno data of samples
    url = 'https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/brca_metabric/data_clinical_sample.txt'
    spl = pd.read_csv(url, sep='\t', index_col=0)
    spl = spl.drop(['#Identifier to uniquely specify a patient.', '#STRING', '#1', 'PATIENT_ID'])
    spl = spl[spl['Sample Type'] == 'Primary']
    # merge pheno data of samples and patients and make sure that there is only one sample per patient
    pheno = pd.merge(pheno, spl, left_index=True, right_index=True)
    assert pd.Series(pheno['Sample Identifier']).is_unique
    assert pd.Series(pheno.index).is_unique
    pheno = pheno.set_index('Sample Identifier')
    # add age_quartile column according to age quartiles in TCGA-BRCA data (values taken from paper)
    pheno.loc[(pheno['Age at Diagnosis'].astype(float) <= 49, 'age_quartile')] = 'lower'
    pheno.loc[(pheno['Age at Diagnosis'].astype(float) > 67, 'age_quartile')] = 'upper'
    # save expression data and phenotype data
    pheno.to_csv(os.path.join(cwd, 'data', 'metabric_brca_pheno.tsv'), sep='\t')
    df.to_csv(os.path.join(cwd, 'data', 'metabric_brca.tsv'), sep='\t')

def download_known_tfs():
    """Saves known human transcription factors obtained from humantfs.ccbr.utoronto.ca in /data.
    """
    cwd = os.getcwd()
    if not os.path.exists(os.path.join(cwd, 'data')):
        os.mkdir(os.path.join(cwd, 'data'))
    df = pd.read_csv('http://humantfs.ccbr.utoronto.ca/download/v_1.01/TFs_Ensembl_v_1.01.txt', delimiter='\t', index_col=0)
    df.to_csv(os.path.join(cwd, 'data', 'regulators.csv'), sep=',')

def download_protein_coding_genes():
    """Saves protein coding genes obtained from http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/locus_groups/protein-coding_gene.txt
    in /data."""
    cwd = os.getcwd()
    if not os.path.exists(os.path.join(cwd, 'data')):
        os.mkdir(os.path.join(cwd, 'data'))
    df = pd.read_csv('http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/locus_groups/protein-coding_gene.txt', delimiter='\t', index_col=False, dtype=str)
    df = df[['ensembl_gene_id', 'symbol']]
    df.to_csv(os.path.join(cwd, 'data', 'protein-coding_gene.csv'), sep=',', index=False)

