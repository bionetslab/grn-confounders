from .NetworkInferenceWrapper import NetworkInferenceWrapper
import sys
import os
import pandas as pd
import numpy as np
import subprocess
import csv
test_suite = os.path.join(os.path.dirname(__file__))
sys.path.append(test_suite)

class GENIE3Wrapper(NetworkInferenceWrapper):

    def _infer_network(self, expression_data, rank):
        """Method to infer a network from expression data using the GENIE3 algorithm.

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
        main = os.getcwd()
        prefix = 'genie3'+str(rank)
        if not os.path.exists(os.path.join(main, 'temp')):
            os.mkdir(os.path.join(main, 'temp'))
        assert os.path.exists(main, 'algorithms', 'GENIE3'), 'download algorithms directory from\
             https://github.com/bionetslab/grn-confounders to use predefined wrappers.'

        # remove columns with zero standard deviation and normalize columns to unit variance
        expression_data = expression_data.loc[:, (expression_data.std() != 0)]
        expression_data = GENIE3Wrapper.normalizeToUnitVariance(expression_data)

        # save expression_data as csv
        expression_data = expression_data.T # R version of genie wants gene x sample data set
        data_path = os.path.join(main, 'temp', f'{prefix}_expression_data.csv')
        expression_data.to_csv(data_path, sep='\t')

        # get regulators and remove such genes that are not present in expression_data
        ktf_path = os.path.join(main, 'data', 'regulators.csv')
        regulators = np.loadtxt(ktf_path, delimiter='\t', dtype=str)
        regulators = regulators[np.isin(regulators, expression_data.index)]
        regulator_path = os.path.join(main, 'temp', f'{prefix}_regulators.csv')
        pd.DataFrame(regulators).to_csv(regulator_path, sep='\t')

        # set output path
        out_path = os.path.join(main, 'temp', f'{prefix}_link_list.csv')

        prog = os.path.join(main, 'algorithms', 'GENIE3', 'GENIE3_R_wrapper.R')
        command = f'Rscript {prog} {prefix}'
        ret = subprocess.run(command, shell=True)

        # get results
        network = pd.read_csv(out_path, sep='\t', index_col=0)
        network = network.rename({'regulatoryGene': 'source', 'targetGene': 'target'}, axis='columns')
        network['type'] = 'directed'
        
        # remove temporary files
        subprocess.call('rm '+str(out_path), shell=True)
        subprocess.call('rm '+str(data_path), shell=True)
        subprocess.call('rm '+str(regulator_path), shell=True)
        
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
            block = self._inferred_networks[i].copy()
            block = block.dropna()
            top_k_edges = []
            state = []
            j = 0
            for j in range(k):
                try:
                    source = block.iloc[j, 0]
                    target = block.iloc[j, 1]
                    state.append('filled')
                    top_k_edges.append((source, target))
                except:
                    state.append('empty'+str(i))
            return set(top_k_edges), state

    @staticmethod
    def normalizeToUnitVariance(df):
        # normalize gene expression data for each gene vector (col) to unit length
        pd.options.mode.chained_assignment = None 
        for col in df:
            if df[col].std() != 0:
                df[col] = df[col]/df[col].std()
            else:
                print('df contains column '+str(col) + ' with std()==0. Normalization would produce nan value. Remove column ' + str(col) + ' before normalization.')
        pd.options.mode.chained_assignment = 'warn'
        return df