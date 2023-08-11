from confinspect import Downloader
import argparse
from confinspect import Selectors

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
                    prog = 'download_TCGA_cohorts',
                    description = 'Download TCGA gene expression data and phenotype data and list of knon transcription factors as list of regulators for GRN methods.')    
    parser.add_argument('-ct', required=False, nargs='*', help='Use TCGA study abbreviations to specify the cohorts to be downloaded. If not set, all cohorts featured in our paper will be downloaded.')
    args = parser.parse_args()

    if not args.ct:
        ct_sels = ['brca_metabric'] + [str(val) for val in list(Selectors.TCGACancerTypeSelector)]
    else:
        ct_sels = list(args.ct)

    for ct_sel in ct_sels:
        print(f'downloading gene expression data and phenotype data for {ct_sel}. This can take a few minutes...')
        if ct_sel == 'brca_metabric':
            Downloader.download_metabric_data()
        else:
            Downloader.download_TCGA_expression_data(ct_sel)
            Downloader.download_TCGA_phenotype_data(ct_sel)

    print('downloading list of known transcription factors...')
    Downloader.download_known_tfs()
    print('downloading list of protein-coding genes...')
    Downloader.download_protein_coding_genes()
