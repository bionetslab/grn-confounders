import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import os
import scipy
import numpy as np
from scipy.stats import chi2_contingency
import glob

def plot_JIs(ct_sel, confs, algs, N_from, N_to, M_from, M_to, k_max):
    cwd = os.getcwd()
    JI_conf = pd.DataFrame(columns = ['confounder', 'cohort', 'algorithm', 'partType', 'partID', 'k', 'mean JI'])
    rnd_counter = 0
    conf_counter = 0
    fro = min([N_from, M_from])
    to = max([N_to, M_to])
    for alg_sel in algs:
        for conf_sel in confs:
            for ct_sel in ct_sels:
                path = glob.glob(os.path.join(cwd, 'results', 'JI')
                for i in range(fro, to):
                    try:
                        filename = 'cb'+'_'+str(i)+'_'+alg_sel+'_'+conf_sel+'_'+ct_sel+'_jaccInd.csv'
                        df = pd.read_csv(os.path.join(path, filename), sep=',', header=0)
                        conf_counter += 1
                        df['partID'] = i
                        df['partType'] = 'conf'
                        df['algorithm'] = alg_sel
                        df['cohort'] = ct_sel
                        df['confounder'] = conf_sel
                        JI_conf = pd.concat([JI_conf, df])
                    except:
                        pass
                    try:
                        filename = 'rnd'+'_'+str(i)+'_'+alg_sel+'_'+conf_sel+'_'+ct_sel+'_jaccInd.csv'
                        df = pd.read_csv(os.path.join(path, filename), sep=',', header=0)
                        rnd_counter += 1
                        df['partID'] = i
                        df['partType'] = 'rnd'
                        df['algorithm'] = alg_sel
                        df['cohort'] = ct_sel
                        df['confounder'] = conf_sel
                        JI_conf = pd.concat([JI_conf, df])
                    except:
                        pass
    print('number of rows corresponding to random partitions:' + str(rnd_counter))
    print('number of rows corresponding to confounder-based partitions:' + str(conf_counter))
    JI_all_type = JI_conf.copy()
    JI_all_type['partType'] = JI_all_type['partType'].replace(['rnd'], 'random')
    JI_all_type['partType'] = JI_all_type['partType'].replace(['conf'], 'confounder-based')
    JI_all_type.rename(columns={'partType': 'partition type', 'algorithm': 'method'}, inplace=True)
    JI_all = JI_all_type[JI_all_type['k']%100 == 10]
    sns.set(font_scale=1.8)
    sns.set_style('whitegrid')
    for cohort in ct_sels:
        JI = JI_all[JI_all['cohort'] == cohort]
        g = sns.FacetGrid(JI, row="method", col="confounder", hue='partition type', 
                        margin_titles=True, xlim=(0,5050), ylim=(0,1), legend_out=False, height=5)
        ax = g.axes[0,0]
        g.map(sns.lineplot, "k", "mean JI",errorbar='sd').add_legend()
        g.fig.suptitle(str(cohort), y=1.04, x=0.66, fontsize='x-large')
        sns.move_legend(ax, "upper left", frameon=False, bbox_to_anchor=(0.0, 1.5))
        g.savefig(os.path.join(cwd, 'res_plots', str(cohort)+'_grid.pdf'))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
                    prog = 'create_JI_lineplots',
                    description = 'Create line plots of the computed mean Jaccard Indices over k.')
    parser.add_argument('-ct_sels', required=True, nargs='+')
    parser.add_argument('-confs', required=True, nargs='+')
    parser.add_argument('-algs', required=True, nargs='+')
    parser.add_argument('-N_from', required=False, type=int, nargs='?', const=0, default=0)
    parser.add_argument('-N_to', required=False, type=int, nargs='?', const=100, default=100)
    parser.add_argument('-M_from', required=False, type=int, nargs='?', const=0, default=0)
    parser.add_argument('-M_to', required=False, type=int, nargs='?', const=10, default=10)
    parser.add_argument('-k_max', required=False, type=int, nargs='?', const=5000, default=5000)
    args = parser.parse_args()
    plot_JIs(args.ct_sel, args.confs, args.algs, args.N_from, args.N_to, args.M_from, args.M_to, args.k_max)