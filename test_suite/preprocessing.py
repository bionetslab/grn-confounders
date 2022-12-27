import numpy as np
import os
import pickle
import pandas as pd
import csv
from sklearn.preprocessing import Normalizer
import Selectors

def normalizeToUnitVariance(df):
    # normalize gene expression data for each gene vector (col) to unit length
    pd.options.mode.chained_assignment = None  # default='warn'
    for col in df:
        if df[col].std() != 0:
            df[col] = df[col]/df[col].std()
        else:
            print('df contains column '+str(col) + ' with std()==0. Normalization would produce nan value. Remove column ' + str(col) + ' before normalization.')
    pd.options.mode.chained_assignment = 'warn'
    return df

def _var_conf_chi(pheno, conf_sel, var, conf_dict):
    conf_partition = Selectors.get_conf_partition(pheno, conf_dict[conf_sel], conf_sel)
    confusion_table = pd.DataFrame()
    for block in conf_partition:
        var_partition = Selectors.get_conf_partition(pheno.loc[block[1]], conf_dict[var], var, min_block_size=0)
        for var_block in var_partition:
            confusion_table.loc[var_block[0], block[0]] = len(var_block[1])
    try:                
        confusion_table = confusion_table.dropna()
        chi2, p, dof, ex = chi2_contingency(confusion_table, correction=False)
    except:
        print('no chi^2 test possible.')
        p = 10000000000
    return p, confusion_table