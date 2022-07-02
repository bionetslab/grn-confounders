import GENIE3Wrapper
import ARACNEWrapper
import NetworkInferenceWrapper
#import Selectors
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), 'GENIE3'))
sys.path.append(os.path.join(os.path.dirname(__file__)))
import preprocessing as prp
import argparse
import numpy as np
import pandas as pd
import csv
import datetime
from sklearn.preprocessing import Normalizer

# partitioning of data set based on confounder
# required: preprocessed data; one col per gene, one row per sample; sample array in same order; gene array in same order
# specify data set path, genes you want to take into account and confounder for partitioning


# only testing
wrapper = ARACNEWrapper.ARACNEWrapper()
#wrapper = GENIE3Wrapper.GENIE3Wrapper()
#df = Selectors.download_TCGA_expression_data(Selectors.CancerTypeSelector.BLCA)

####################################################################
# parameters
####################################################################

cancer_type = 'TCGA-BLCA'
confounder = 'race'
do_preprocessing = False

####################################################################
# path handling
####################################################################
dateStart = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
print('Result folder time stamp: ' + str(dateStart))
cwd = os.path.join(os.path.dirname(__file__))
raw_data_dir = os.path.join(cwd, '..', 'input', str(cancer_type)+'_raw_data')
data_dir = os.path.join(cwd, '..', 'input', str(cancer_type)+'_data')
results_dir = os.path.join(cwd, 'results')
timestamp_dir = os.path.join(results_dir, str(dateStart))
conf_results_dir = os.path.join(timestamp_dir, str(cancer_type)+'_conf_results')
rnd_results_dir = os.path.join(timestamp_dir, str(cancer_type)+'_rnd_results')

if not os.path.exists(results_dir):
    os.mkdir(results_dir)
if not os.path.exists(timestamp_dir):
    os.mkdir(timestamp_dir)
if not os.path.exists(conf_results_dir):
    os.mkdir(conf_results_dir)
if not os.path.exists(rnd_results_dir):
    os.mkdir(rnd_results_dir)

####################################################################
# start of script
####################################################################

if do_preprocessing:
    prp.preprocessToCsv(raw_data_dir, data_dir, cancer_type) # cancer_type only needed for file naming; TODO: later also for data download

####################################################################
# prepare final data set
####################################################################

######################
# get data: one col per gene, one row per sample and remove columns col where col.std() is 0
######################
df = pd.read_csv(os.path.join(data_dir, str(cancer_type)+'_data_set.csv'), sep='\t', index_col=0, header=0)
df = df.iloc[:, :203] # TODO remove later
df = df.loc[:, (df.std() != 0)]
for col in df:
    df.rename({col: col.split('.')[0]}, axis=1, inplace=True)

prp.testDuplicateGenes(df.columns.copy())

######################
# filter samples (rows): only keep 'Primary Tumor' samples
######################
tum_samples = prp.getPrimaryTumorIndices(os.path.join(raw_data_dir, str(cancer_type)+'.GDC_phenotype.tsv'))[1]
df = df.filter(items = tum_samples, axis=0)

######################
# get list of regulators
######################

#regulators = np.genfromtxt(fname=os.path.join(raw_data_dir, 'TFs_Ensembl_v_1.01.txt'), delimiter="\t", dtype=str)

######################
# group sample ids based on confounder-induced partitions
######################
(blocknames, conf_partition) = prp.getSampleIDsByConfounder(os.path.join(raw_data_dir, str(cancer_type)+'.GDC_phenotype.tsv'), confounder)

reg = pd.read_csv('http://humantfs.ccbr.utoronto.ca/download/v_1.01/TFs_Ensembl_v_1.01.txt', delimiter='\t', index_col=0)
reg.to_csv(os.path.join(cwd, '..', 'data', 'regulators.csv'), sep='\t')

wrapper.partition = conf_partition
#print(df.head())
wrapper.expression_data = df
print(df.head())

wrapper.infer_networks()
print(wrapper._inferred_networks[0].head())
s = wrapper._get_top_k_edges(0, 5)
print(s)
