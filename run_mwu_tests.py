import argparse
from scipy.stats import mannwhitneyu
from test_suite import Selectors
from test_suite import preprocessing
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns

def run_mwu_tests(ct_sel, conf_sels, alg_sels, k, N_from, N_to, M_from, M_to):
    """For each k in range(10, k, 100), compute the one-sided Mann-Whitney-U test of the random-based distribution and cnfounder-based
    distribution of mean Jaccard Indices under the alternative hypothesis of smaller values in the confounder-based distribution.
    Parameters
    ----------
    ct_sel : str
        Cohort to test.
    conf_sels : list
        Confounders underlying the networks to be tested.
    alg_sels : list
        Algorithms that produced the results to be tested.
    k : int
        Mean Jaccard Index ill be computed for range(10, k, 100).
    N_from : int
        Starting index of random-based networks to be tested. Must be equal in G_all, random, and confounder-based networks.
    N_to : int
        Stopping index of random-based networks to be tested. Must be equal in G_all, random, and confounder-based networks.
    M_from : int
        Starting index of confounder-based networks to be tested. Must be equal in G_all, random, and confounder-based networks.
    M_to : int
        Stopping index of confounder-based networks to be tested. Must be equal in G_all, random, and confounder-based networks.
    """
    mwu_fractions = pd.DataFrame(columns=['cohort', 'confounder/variable', 'method', 'frac'])
    for conf_sel in conf_sels:
        for alg_sel in alg_sels:
            mwu_fractions = pd.concat([mwu_fractions, run_mwu_tests_single_test(ct_sel, conf_sel, alg_sel, k, N_from, N_to, M_from, M_to)])         
    mwu_fractions.to_csv(os.path.join('results', ct_sel + '_mwu_fractions.csv'))
    visualize_mwus(ct_sel, mwu_fractions, save=True)

def run_mwu_tests_single_test(ct_sel, conf_sel, alg_sel, k, N_from, N_to, M_from, M_to):
    """For each k in range(10, k, 100), compute the one-sided Mann-Whitney-U test of the random-based distribution and cnfounder-based
    distribution of mean Jaccard Indices under the alternative hypothesis of smaller values in the confounder-based distribution.
    Visualize results in sarm plot, grouped by method and confounder. Run tests for a single confounder and algorithm.
    Parameters
    ----------
    ct_sel : str
        Cohort to test.
    conf_sel : str
        Confounder underlying the networks to be tested.
    alg_sel : str
        Algorithm that produced the results to be tested.
    k : int
        Mean Jaccard Index ill be computed for range(10, k, 100).
    N_from : int
        Starting index of random-based networks to be tested. Must be equal in G_all, random, and confounder-based networks.
    N_to : int
        Stopping index of random-based networks to be tested. Must be equal in G_all, random, and confounder-based networks.
    M_from : int
        Starting index of confounder-based networks to be tested. Must be equal in G_all, random, and confounder-based networks.
    M_to : int
        Stopping index of confounder-based networks to be tested. Must be equal in G_all, random, and confounder-based networks.
    """
    rep_k = range(10, k, 100)
    mwu_fractions = pd.DataFrame(columns=['cohort', 'confounder/variable', 'method', 'frac'])
    path = os.path.join(os.getcwd(), 'results', 'JI')
    c = pd.DataFrame(columns=['size intersection', 'size union', 'state', 'k','mean JI'])
    r = pd.DataFrame(columns=['size intersection', 'size union', 'state', 'k','mean JI'])
    for i in range(M_from, M_to):
        filename = 'cb'+'_'+str(i)+'_'+alg_sel+'_'+conf_sel+'_'+ct_sel+'_jaccInd.csv'
        assert os.path.exists(os.path.join(path, filename)), 'JIs from repetition ' + str(i) + ' on confounder-based partition do not exist. Run run_tests.py first, and include the missing repetition in the range M_from, M_to.'
        c = pd.concat([c, pd.read_csv(os.path.join(path, filename), header=0)])
    for i in range(N_from, N_to):
        filename = 'rnd'+'_'+str(i)+'_'+alg_sel+'_'+conf_sel+'_'+ct_sel+'_jaccInd.csv'
        assert os.path.exists(os.path.join(path, filename)), 'JIs from random partition ' + str(i) + ' do not exist. Run run_tests.py first, and include the missing partition in the range N_from, N_to.'
        r = pd.concat([r, pd.read_csv(os.path.join(path, filename), header=0)])
    significants = 0
    for k in rep_k:
        r_k = r[r['k'] == k]
        c_k = c[c['k'] == k]
        mwu = mannwhitneyu(c_k['mean JI'].astype(float),r_k['mean JI'].astype(float), alternative='less')
        if mwu.pvalue < 0.05:
            significants += 1
    mwu_fractions = pd.concat([mwu_fractions, pd.DataFrame({'cohort':[ct_sel], 'confounder/variable':[str(conf_sel)], 
                        'method':[str(alg_sel)], 'frac':[significants/len(rep_k)]})])
    return mwu_fractions
    
def visualize_mwus(ct_sel, mwu_fractions, save=True):
    """Visualize fraction of significant k in sarm plot, grouped by method and confounder.
    Parameters
    ----------
    ct_sel : str
        String identifier of cancer type to be tested.
    mwu_fractions : pd.DataFrame
        pd.DataFrame containing fractions of significant k per method and per confounder.
    save : bool
        Internal boolean whether to save the lineplots or not. Defaults to True.
    """
    plt.tight_layout()
    sns.set(font_scale=1.2)
    sns.set_style("darkgrid")
    h, axs_z = plt.subplots()
    h.tight_layout(pad=5.0)

    cur = mwu_fractions[mwu_fractions['cohort'] == ct_sel]
    cur_var = pd.DataFrame({'method':cur['confounder/variable'], 'confounder/variable':cur['method'], 'frac':cur['frac']})
    cur = pd.concat([cur, cur_var])

    g = sns.swarmplot(ax=axs_z, x="method", y="frac",hue="confounder/variable",data=cur, legend=True, s=15)
    g.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), ncol=1)
    ax = g.axes
    ax.yaxis.grid(True)
    ax.xaxis.grid(True)

    ax.set_title(ct_sel, loc='left', fontdict={'fontsize': 28,'verticalalignment': 'baseline','horizontalalignment': 'left'})
    g.set(xlabel='', ylabel = "fraction of k with \n significant p-value",title=ct_sel)
    ax.yaxis.label.set_size(18)
    g.tick_params(labelsize=16)
    plt.yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])

    if save:
        h.savefig(os.path.join(cwd, 'mwu_comparison', str(ct_sel)+'_mwu_comparison.pdf'))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
                    prog = 'run_mwu_tests',
                    description = 'Run Mann-Whitney-U Test to compare random-based and confounder-based distributions of mean Jaccard Indices at k, for all k. Obtain fraction of k with significant p-value.')
    parser.add_argument('-ct', required=True)
    parser.add_argument('-conf', required=True, nargs='+')
    parser.add_argument('-alg', required=True, nargs='+', choices=[str(sel) for sel in list(Selectors.AlgorithmSelector)])
    parser.add_argument('-N_from', required=False, type=int, nargs='?', const=0, default=0)
    parser.add_argument('-M_from', required=False, type=int, nargs='?', const=0, default=0)
    parser.add_argument('-N_to', required=False, type=int, nargs='?', const=20, default=20)
    parser.add_argument('-M_to', required=False, type=int, nargs='?', const=20, default=20)
    parser.add_argument('-k', required=False, type=int, nargs='?', const=5000, default=5000)
    args = parser.parse_args()
    assert args.N_from < args.N_to
    assert args.M_from < args.M_to
    run_mwu_tests(args.ct, args.conf, args.alg, args.k, args.N_from, args.N_to, args.M_from, args.M_to)
