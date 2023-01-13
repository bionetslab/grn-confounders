from scipy.stats import wilcoxon
from test_suite import Selectors
from test_suite import preprocessing
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

def compute_and_plot_JIs(ct_sel, conf_sels, alg_sels, fro, to, k, block_ids):
    print('Compute and plot JIs of G_all with random and confounder-linked networks of each block.')
    rep_k = range(10, k, 100)
    for conf_sel in conf_sels:
        for alg_sel in alg_sels:
            if conf_sel != Selectors.ConfounderSelector.NONE:
                compute_and_plot_JIs_single(ct_sel, conf_sel, alg_sel, fro, to, rep_k, block_ids)

def compute_and_plot_JIs_single(ct_sel, conf_sel, alg_sel, fro, to, rep_k, block_ids):
    algorithm_wrapper = Selectors.get_algorithm_wrapper(Selectors.AlgorithmSelector(alg_sel))
    print('Compute mean JIs for G_all and confounder-based networks...')
    networks_blocks_conf = {bl: pd.DataFrame(columns=['size intersection', 'size union', 'state', 'k', 'mean JI', 'partID']) for bl in block_ids}
    conf_results = {l:{bl:[] for bl in block_ids} for l in range(fro, to)}
    for j in range(fro, to):
        plain = pd.read_csv(os.path.join(os.getcwd(), 'results', 'networks', 'conf_part'+str(j)+'_all_'+alg_sel+'_'+ct_sel+'_NONE_gene_list.csv'), header=0)
        for bl in block_ids:
            block_network = os.path.join(os.getcwd(), 'results', 'networks', 'conf_part'+str(j)+'_'+bl+'_'+alg_sel+'_'+ct_sel+'_'+conf_sel+'_gene_list.csv')
            algorithm_wrapper._inferred_networks = [plain, pd.read_csv(block_network, header=0)]
            for k in rep_k:
                ji, state, s_int, s_un = algorithm_wrapper.mean_jaccard_index_at_k(k)
                conf_results[j][bl].append(ji)
            JI_conf = pd.DataFrame({'k': rep_k, 'mean JI': conf_results[j][bl], 'partID':j,'confounder':conf_sel, 'method':alg_sel, 'cohort':ct_sel})
            networks_blocks_conf[bl] = pd.concat([networks_blocks_conf[bl], JI_conf])

    print('Compute mean JIs for G_all and random networks...')
    networks_blocks_rnd = {bl: pd.DataFrame(columns=['size intersection', 'size union', 'state', 'k', 'mean JI', 'partID']) for bl in block_ids}
    rnd_results = {l:{bl:[] for bl in block_ids} for l in range(fro, to)}
    for j in range(fro, to):
        plain = pd.read_csv(os.path.join(os.getcwd(), 'results', 'networks', 'conf_part'+str(j)+'_all_'+alg_sel+'_'+ct_sel+'_NONE_gene_list.csv'), header=0)
        for bl in block_ids:
            block_network = os.path.join(os.getcwd(), 'results', 'networks', 'rnd_part'+str(j)+'_'+bl+'_'+alg_sel+'_'+ct_sel+'_'+conf_sel+'_gene_list.csv')
            algorithm_wrapper._inferred_networks = [plain, pd.read_csv(block_network, header=0)]
            for k in rep_k:
                ji, state, s_int, s_un = algorithm_wrapper.mean_jaccard_index_at_k(k)
                rnd_results[j][bl].append(ji)
            JI_rnd = pd.DataFrame({'k': rep_k, 'mean JI': rnd_results[j][bl], 'partID':j, 'confounder':conf_sel,
                                    'method':alg_sel, 'cohort':ct_sel})
            networks_blocks_rnd[bl] = pd.concat([networks_blocks_rnd[bl], JI_rnd])

    print('Plot JIs...')
    for bl in block_ids:
        visualize_JIs(networks_blocks_conf[bl], networks_blocks_rnd[bl], bl)

def visualize_JIs(block_conf, block_rnd, bl, save=False):
    plt.figure()
    sns.set(font_scale=1.2)
    sns.set_style('whitegrid')
    block_conf['partition type'] = str(bl)
    block_rnd['partition type'] = 'random with size\nof ' + str(bl) + ' block'
    JI = pd.concat([block_conf, block_rnd])
    g = sns.lineplot(JI, x='k', y='mean JI', hue='partition type', errorbar='sd')
    g.set_title(alg_sel + ' ' + ct_sel + ' ' + conf_sel)
    plt.legend(bbox_to_anchor=(1.02, 0.9), loc='upper left', title='block compared with $G^{all}$', borderaxespad=0)
    if save:
        g.figure.savefig(os.path.join(os.getcwd(), 'plots', alg_sel+'_'+conf_sel+'_'+ct_sel+'_compare_conf_and_rnd.pdf'))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
                    prog = 'compare_g_all',
                    description = 'Compare networks inferred from blocks with G_all, then compare netorks from random blocks with G_all, by computing mean JIs. Then plot the JIs.')
    parser.add_argument('-ct', required=True)
    parser.add_argument('-conf', required=True, nargs='+')
    parser.add_argument('-alg', required=True, nargs='+', choices=[str(sel) for sel in list(Selectors.AlgorithmSelector)])
    parser.add_argument('-fro', required=False, type=int, nargs='?', const=0, default=0)
    parser.add_argument('-to', required=False, type=int, nargs='?', const=20, default=20)
    parser.add_argument('-k', required=False, type=int, nargs='?', const=5000, default=5000)
    parser.add_argument('-block_ids', required=True, nargs='+')
    args = parser.parse_args()
    assert args.fro < args.to
    compute_and_plot_JIs(args.ct, args.conf, args.alg, args.fro, args.to, args.k, args.block_ids)
    
