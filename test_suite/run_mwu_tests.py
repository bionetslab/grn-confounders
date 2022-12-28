from scipy.stats import mannwhitneyu
import Selectors
import preprocessing
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns

def run_mwu_tests(ct_sel, conf_sels, alg_sels, rep_k, N_from, N_to, M_from, M_to):
    mwu_fractions = pd.DataFrame(columns=['cohort', 'confounder/variable', 'method', 'frac'])
    for conf_sel in conf_sels:
        for alg_sl in alg_sels:
            mwu_fractions = pd.concat(mwu_fractions, run_mwu_tests_single_test(ct_sel, conf_sel, alg_sel, rep_k, N_from, N_to, M_from, M_to))         
    mwu_fractions.to_csv(os.path.join('results', ct_sel + '_mwu_fractions.csv'))
    visualize_mwus(ct_sel, mwu_fractions)

def run_mwu_tests_single_test(ct_sel, conf_sel, alg_sel, rep_k, N_from, N_to, M_from, M_to):
    mwu_fractions = pd.DataFrame(columns=['cohort', 'confounder/variable', 'method', 'frac'])
    path = os.path.join(os.getcwd(), 'results', 'JI')
    c = pd.DataFrame(columns=['size intersection', 'size union', 'state', 'k','mean JI'])
    r = pd.DataFrame(columns=['size intersection', 'size union', 'state', 'k','mean JI'])
    for i in range(M_from, M_to):
        filename = 'cb'+'_'+str(i)+'_'+alg_sel+'_'+conf_sel+'_'+ct_sel+'_jaccInd.csv'
        c = pd.concat(c, pd.read_csv(os.path.join(path, filename), sep=',', header=0))
    for i in range(N_from, N_to):
        filename = 'rnd'+'_'+str(i)+'_'+alg_sel+'_'+conf_sel+'_'+ct_sel+'_jaccInd.csv'
        r = pd.concat(r, pd.read_csv(os.path.join(path, filename), sep=',', header=0))
    significants = 0
    for k in rep_k:
        r_k = r[r['k'] == k]
        c_k = c[c['k'] == k]
        mwu = mannwhitneyu(c_k['mean JI'],r_k['mean JI'], alternative='less')
        if mwu.pvalue < 0.05:
            significants += 1
    mwu_fractions = pd.concat([mwu_fractions, pd.DataFrame({'cohort':ct_sel, 'confounder/variable':str(conf_sel), 
                        'method':str(alg_sel), 'frac':significants/len(rep_k)})])
    return mwu_fractions

    def visualize_mwus(ct_sel, mwu_fractions, save=False):
        plt.tight_layout()
        sns.set(font_scale=1.2)
        sns.set_style("darkgrid")
        palette={'age':"#2a60e1", 'ethnicity':"#e8e049", 'sex':"#69d447", 'WGCNA':'#15a5d9', 'CEMiTool':'#f9a62d', 'ARACNe-AP':'#d54f13'}
        h, axs_z = plt.subplots()
        h.tight_layout(pad=5.0)

        cur = mwu_fractions[mwu_fractions['cohort'] == ct_sel]
        cur_var = pd.DataFrame({'method':cur['confounder/variable'], 'confounder/variable':cur['method'], 'frac':cur['frac']})
        cur = pd.concat([cur, cur_var])

        g = sns.swarmplot(ax=axs_z, x="method", y="frac",hue="confounder/variable",data=_cur, legend=True, s=15, palette=palette)
        g.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), ncol=1)
        ax = g.axes
        ax.yaxis.grid(True)
        ax.xaxis.grid(True)
        ax.axvline(2.5, ls='--')

        ax.set_title(ct_sel, loc='left', fontdict={'fontsize': 28,'verticalalignment': 'baseline','horizontalalignment': 'left'})
        g.set(xlabel='', ylabel = "fraction of k with \n significant p-value",title=ct_sel)
        ax.yaxis.label.set_size(18)
        g.tick_params(labelsize=16)
        plt.yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])

        if save:
            h.savefig(os.path.join(cwd, 'mwu_comparison', str(ct_sel)+'_mwu_comparison.pdf'))

if __name__ == '__main__':
    # TODO parser
    run_mwu_tests(ct_sel, conf_sels, alg_sels, rep_k, N_from, N_to, M_from, M_to)
