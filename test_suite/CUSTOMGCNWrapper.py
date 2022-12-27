from NetworkInferenceWrapper import NetworkInferenceWrapper
import sys
import os
import pandas as pd
import numpy as np
import subprocess
from . import preprocessing as prp
import csv

test_suite = os.path.join(os.path.dirname(__file__))
sys.path.append(test_suite)


class CUSTOMGCNWrapper(NetworkInferenceWrapper):

    def _infer_network(self, expression_data, rank):
        """Method to infer a network from expression data using a custom GCN inference method.

        Parameters
        ----------
        expression_data : pd.DataFrame
            Gene expression data stored in a data frame with sample identifiers as indices and
            gene symbols as column names.

        Returns
        -------
        inferred_network : pd.DataFrame
            A data frame whose entries in the column 'score' correspond to edge scores in inferred network. For undirected networks,
            columns 'node_lower' and 'node_upper' contain the gene symbols of the nodes that are connected by the edge.
            Fr directed networks, these columns are named 'source' and 'target'.
        """
        main = os.path.join(test_suite, '..')
        prefix = 'customgcn'+str(rank)

        data_path = os.path.join(main, 'temp', f'{prefix}_expression_data.csv')

        # add preprocessing, if required

        # save expression data
        expression_data.to_csv(data_path, sep='\t')

        out_path = os.path.join(main, 'temp', f'{prefix}_edge_list.csv')

        command = f'call_method_here {prefix}'
        subprocess.run(command, shell=True)

        # get results
        network = pd.read_csv(out_path, sep='\t', index_col=0)
        network = network.sort_values(by=['score', 'source', 'target'], axis=0, ascending=False)
        network['type'] = 'undirected'

        # remove temporary files
        subprocess.call('rm ' + str(out_path), shell=True)
        subprocess.call('rm ' + str(data_path), shell=True)

        return network

    def _get_top_k_edges(self, i, k):
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
        block = self._inferred_networks[i]
        top_k_edges = []
        state = []
        for j in range(k):
            try:
                gene_1 = block.iloc[j, 0]
                gene_2 = block.iloc[j, 1]
                state.append('filled')
                if(gene_1 < gene_2):
                    top_k_edges.append((gene_1, gene_2))
                else:
                    top_k_edges.append((gene_2, gene_1))
            except:
                state.append('empty' + str(i))
        return set(top_k_edges), state
