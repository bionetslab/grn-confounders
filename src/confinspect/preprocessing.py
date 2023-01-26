import os
import pandas as pd

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
