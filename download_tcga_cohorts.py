from confinspect import Downloader
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
                    prog = 'download_TCGA_cohorts',
                    description = 'Download TCGA gene expression data and phenotype data and list of knon transcription factors as list of regulators for GRN methods.')    
    parser.add_argument('-tfs', action='store_true', help='If set, download list of human known transcription factors')
    parser.add_argument('-ct', required=True, nargs='+', help='Use TCGA study abbreviations to specify the cohorts to be downloaded')
    args = parser.parse_args()
    for ct_sel in list(args.ct):
        print(f'downloading gene expression data for {ct_sel}. This can take a few minutes...')
        Downloader.download_TCGA_expression_data(ct_sel)
        print(f'downloading pheno type data for {ct_sel}. This can take a few minutes...')
        Downloader.download_TCGA_phenotype_data(ct_sel)
    if args.tfs:
        print('downloading list of known transcription factors...')
        Downloader.download_known_tfs()
