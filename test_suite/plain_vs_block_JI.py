from abc import ABC, abstractmethod
import itertools as itt
import os
import pandas as pd
from glob import glob

def mean_jaccard_index_at_k(k, inferred_networks):
    """Returns the mean Jaccard index for the top k edges in the inferred networks.

    Parameters
    ----------
    k : int
        Number of edges which should be considered when computing Jaccard indices.

    Returns
    -------
    mean_jaccard_index : float
        The mean Jaccard index across of pairwise comparisons of the top k edges in the 
        networks inferred for the blocks of the partition.

    state_k: str
        indicating which networks contained at least k edges.

    size_intersection: int
        Size of the intersection at k.

    size_union: int
        Size of the union at k.
    """
    sum_jaccard_indices = 0.0
    num_comparisons = 0
    state_k = ''
    for i, j in itt.combinations(range(len(inferred_networks)), 2):
        top_k_edges_i, state_i = get_top_k_edges(i, k, inferred_networks)
        top_k_edges_j, state_j = get_top_k_edges(j, k, inferred_networks)
        if 'empty'+str(i) in state_i:
            state_k += 'empty'+str(i) + '_'
        if 'empty'+str(j) in state_j:
            state_k += 'empty'+str(j) + '_'
        if not ('empty'+str(i) in state_i or 'empty'+str(j) in state_j):
            state_k += 'filled'+str(i)+str(j)+'_'

        size_intersection = len(top_k_edges_i.intersection(top_k_edges_j))
        size_union = len(top_k_edges_i.union(top_k_edges_j))
        if size_union > 0:
            sum_jaccard_indices += size_intersection / size_union
        else:
            sum_jaccard_indices += 0
        num_comparisons += 1
    return sum_jaccard_indices / num_comparisons, state_k, size_intersection, size_union

def get_top_k_edges(i, k, inferred_networks, alg_sel):
    """Method to return the top k edges for the inferred network for block i.

    Parameters
    ----------
    i : int
        ID of partion block. Must fall in range(0, len(self.partition)).
    k : Maximal number of edges to be returned.

    Returns
    -------
    top_k_edges : set
        Set of tuples encoding edges. For edges without a sense, use tuples of form (<gene_1>, <gene_2>),
        where <gene_1> and <gene_2> are gene symbols. For edges with a sense (e.g., positive or negative
        correlation), use tuples of form (<gene_1>, <gene_2>, <sense>), where <sense> is either -1 or 1.
        For undirected edges, ensure that <gene_1> <= <gene_2> for all tuples contained in edge set.
    """
    block = inferred_networks[i]
    top_k_edges = []
    state = []
    for j in range(k):
        try:
            gene_1 = block.iloc[j, 0]
            gene_2 = block.iloc[j, 1]
            state.append('filled')
            if alg_sel == 'CEMI' or alg_sel == 'WGCNA':
                if(gene_1 < gene_2):
                    top_k_edges.append((gene_1, gene_2))
                else:
                    top_k_edges.append((gene_2, gene_1))
            elif alg_sel == 'ARACNE':
                top_k_edges.append((source, target))
        except:
            state.append('empty' + str(i))
    return set(top_k_edges), state

if __name__ == "__main__":
    alg_sels = ['WGCNA', 'CEMI']
    ct_sels = ['KIRC', 'LUSC', 'KIRP']
    conf_sels = ['sex']
    for alg_sel in alg_sels:
        for ct_sel in ct_sels:
            for conf_sel in conf_sels:
                conf_results = {ct_sel:{conf_sel:{alg_sel:{l:{bl:[] for bl in ['0', '1']} for l in range(0,10)}}}}
                # here 1 if WGCNA or CEMI on conf vs entire, 20 else
                for j in range(10, 20):
                    print('testing blocks in partition '+str(j))
                    # here j = 0 on none-entire-network if WGCNA, CEMI
                    if alg_sel == 'WGCNA' or alg_sel == 'CEMI':
                        plain = pd.read_csv(os.path.join(os.getcwd(), alg_sel+'_RUN', 'results', 'networks', 'conf_part'+str(0)+'_block0_'+alg_sel+'_'+ct_sel+'_none_gene_list.csv'), header=0)
                    elif alg_sel == 'ARACNE':
                        plain = pd.read_csv(os.path.join(os.getcwd(), alg_sel+'_RUN', 'results', 'networks', 'conf_part'+str(j)+'_block0_'+alg_sel+'_'+ct_sel+'_none_gene_list.csv'), header=0)
                    for bl in ['0', '1']:
                        # here conf if conf vs entire, here rnd if rnd block vs entire
                        block0 = os.path.join(os.getcwd(), alg_sel+'_RUN', 'results', 'networks', 'rnd_part'+str(j)+'_block'+bl+'_'+alg_sel+'_'+ct_sel+'_'+conf_sel+'_gene_list.csv')
                        inferred_networks = [plain, pd.read_csv(block0, header=0)]
                        index=[]
                        network_state = []
                        intersections = []
                        unions = []
                        for k in range(10, 5000, 50):
                            ji, state, s_int, s_un = mean_jaccard_index_at_k(k, inferred_networks, alg_sel)
                            conf_results[ct_sel][conf_sel][alg_sel][j][bl].append(ji)
                            #index.append(k)
                            unions.append(s_un)
                            intersections.append(s_int)
                            network_state.append(state)
                        pd.DataFrame({'size intersection': intersections, 'size union': unions, 'state': network_state, 'k': range(10, 5000, 50), 'mean JI': conf_results[ct_sel][conf_sel][alg_sel][j][bl]}).to_csv(os.path.join(os.getcwd(), 'results', f'rnd_{j}_{str(alg_sel)}_{str(conf_sel)}_{str(ct_sel)}_{bl}_jaccInd.csv'), index=False)
