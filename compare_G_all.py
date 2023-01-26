import confinspect
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

def run_G_all_tests(ct_sel, conf_sels, alg_sels, fro, to, k, block_ids):
    """Compare G_all with networks inferred from single confounder-based and random blocks. Networks must be saved in results/networks/
    as conf_partj_all_alg_sel_ct_sel_NONE_gene_list.csv (G_all) and conf_partj_all_alg_sel_ct_sel_conf_sel_gene_list.csv
    (confounder based), and rnd_partj_all_alg_sel_ct_sel_conf_sel_gene_list.csv (random based), for j=[fro,...,to).
    Parameters
    ----------
    ct_sel : str
        Cohort to test.
    conf_sels : list
        Confounders underlying the networks to be tested.
    alg_sels : list
        Algorithms that produced the results to be tested.
    fro : int
        Starting index of networks to be tested. Must be equal in G_all, random, and confounder-based networks.
    to : int
        Stopping index of networks to be tested. Must be equal in G_all, random, and confounder-based networks.
    k : int
        Mean Jaccard Index ill be computed for range(10, k, 100).
    block_ids : list
        String identifiers of the confounder-based blocks underlying the networks to be tested. Order is irrelevant.
    """
    for conf_sel in conf_sels:
        for alg_sel in alg_sels:
            if conf_sel != confinspect.Selectors.ConfounderSelector.NONE:
                run_G_all_test_single(ct_sel, conf_sel, alg_sel, fro, to, k, block_ids)

def run_G_all_test_single(ct_sel, conf_sel, alg_sel, fro, to, k, block_ids):
    """Compare G_all with networks inferred from single confounder-based and random blocks. Networks must be saved in results/networks/
    as conf_partj_all_alg_sel_ct_sel_NONE_gene_list.csv (G_all) and conf_partj_all_alg_sel_ct_sel_conf_sel_gene_list.csv
    (confounder based), and rnd_partj_all_alg_sel_ct_sel_conf_sel_gene_list.csv (random based), for j=[fro,...,to). Run tests for single
    confounder and algorithm.
    Parameters
    ----------
    ct_sel : str
        Cohort to test.
    conf_sel : str
        Confounder underlying the networks to be tested.
    alg_sel : str
        Algorithm that produced the results to be tested.
    fro : int
        Starting index of networks to be tested. Must be equal in G_all, random, and confounder-based networks.
    to : int
        Stopping index of networks to be tested. Must be equal in G_all, random, and confounder-based networks.
    k : int
        Mean Jaccard Index ill be computed for range(10, k, 100).
    block_ids : list
        String identifiers of the confounder-based blocks underlying the networks to be tested. Order is irrelevant.
    """
    rep_k = range(10, k, 100)
    algorithm_wrapper = confinspect.Selectors.get_algorithm_wrapper(confinspect.Selectors.AlgorithmSelector(alg_sel))
    print('Compute mean JIs for G_all and confounder-based networks, and for G_all and random-based networks...')
    networks_blocks_conf = {bl: pd.DataFrame(columns=['size intersection', 'size union', 'state', 'k', 'mean JI', 'partID']) for bl in block_ids}
    conf_results = {l:{bl:[] for bl in block_ids} for l in range(fro, to)}
    networks_blocks_rnd = {bl: pd.DataFrame(columns=['size intersection', 'size union', 'state', 'k', 'mean JI', 'partID']) for bl in block_ids}
    rnd_results = {l:{bl:[] for bl in block_ids} for l in range(fro, to)}
    for j in range(fro, to):
        plain = pd.read_csv(os.path.join(os.getcwd(), 'results', 'networks', 'conf_part'+str(j)+'_all_'+alg_sel+'_'+ct_sel+'_NONE_gene_list.csv'), header=0)
        for bl in block_ids:
            for k in rep_k:
                appendix = str(j)+'_'+bl+'_'+alg_sel+'_'+ct_sel+'_'+conf_sel+'_gene_list.csv'
                # compare confounder-based network to G_all
                block_network = os.path.join(os.getcwd(), 'results', 'networks', 'conf_part' + appendix)
                algorithm_wrapper._inferred_networks = [plain, pd.read_csv(block_network, header=0)]
                ji, state, s_int, s_un = algorithm_wrapper.mean_jaccard_index_at_k(k)
                conf_results[j][bl].append(ji)
                # compare randombased network to G_all
                block_network = os.path.join(os.getcwd(), 'results', 'networks', 'rnd_part' + appendix)
                algorithm_wrapper._inferred_networks = [plain, pd.read_csv(block_network, header=0)]
                ji, state, s_int, s_un = algorithm_wrapper.mean_jaccard_index_at_k(k)
                rnd_results[j][bl].append(ji)
            JI_conf = pd.DataFrame({'k': rep_k, 'mean JI': conf_results[j][bl], 'partID':j,'confounder':conf_sel, 
                                    'method':alg_sel, 'cohort':ct_sel})
            networks_blocks_conf[bl] = pd.concat([networks_blocks_conf[bl], JI_conf])
            JI_rnd = pd.DataFrame({'k': rep_k, 'mean JI': rnd_results[j][bl], 'partID':j, 'confounder':conf_sel,
                                    'method':alg_sel, 'cohort':ct_sel})
            networks_blocks_rnd[bl] = pd.concat([networks_blocks_rnd[bl], JI_rnd])
    visualize_JIs(networks_blocks_conf, networks_blocks_rnd)

def visualize_JIs(block_conf, block_rnd, save=True):
    """Visualize comparison of G_all with networks inferred from single confounder-based and random blocks in line plots over k.
    Parameters
    ----------
    block_conf : dict
        Dict containing DataFrame per block with mean Jaccard Indices over k and information about confounder-based blocks: 'k', 'mean JI', 
        'partID', 'confounder', 'method', and 'cohort'.
    block_rnd : dict
        Dict containing DataFrame with mean Jaccard Indices over k and information about random blocks: 'k', 'mean JI', 'partID',
        'confounder', 'method', and 'cohort'.
    save : bool
        Internal boolean whether to save the lineplots or not. Defaults to True.
    """
    plt.figure()
    sns.set(font_scale=1.2)
    sns.set_style('whitegrid')
    for bl in block_ids:
        c = block_conf[bl]
        r = block_rnd[bl]
        c['partition type'] = str(bl)
        r['partition type'] = 'random with size\nof ' + str(bl) + ' block'
        JI = pd.concat([c, r])
        g = sns.lineplot(JI, x='k', y='mean JI', hue='partition type', errorbar='sd')
        g.set_title(alg_sel + ' ' + ct_sel + ' ' + conf_sel)
    plt.legend(bbox_to_anchor=(1.02, 0.9), loc='upper left', title='block underlying network compared with $G^{all}$', borderaxespad=0)
    if save:
        g.figure.savefig(os.path.join(os.getcwd(), 'plots', alg_sel+'_'+conf_sel+'_'+ct_sel+'_compare_conf_and_rnd.pdf'))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
                    prog = 'compare_g_all',
                    description = 'Compare networks inferred from blocks with G_all, then compare netorks from random blocks with G_all. Plot mean JIs in lineplot over k.')
    parser.add_argument('-ct', required=True)
    parser.add_argument('-conf', required=True, nargs='+')
    parser.add_argument('-alg', required=True, nargs='+', choices=[str(sel) for sel in list(confinspect.Selectors.AlgorithmSelector)])
    parser.add_argument('-fro', required=False, type=int, nargs='?', const=0, default=0)
    parser.add_argument('-to', required=False, type=int, nargs='?', const=20, default=20)
    parser.add_argument('-k', required=False, type=int, nargs='?', const=5000, default=5000)
    parser.add_argument('-block_ids', required=True, nargs='+')
    args = parser.parse_args()
    assert args.fro < args.to
    run_G_all_tests(args.ct, args.conf, args.alg, args.fro, args.to, args.k, args.block_ids)
    
