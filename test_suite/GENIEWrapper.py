from NetworkInferenceWrapper import NetworkInferenceWrapper
import sys
import os
import pandas as pd
import numpy as np
import subprocess
import preprocessing as prp
import csv
test_suite = os.path.join(os.path.dirname(__file__))
sys.path.append(test_suite)

# TODO: define path name conventions

class GENIEWrapper(NetworkInferenceWrapper):

    def _infer_network(self, expression_data):

        """Method to infer a network from expression data using the GENIE3 algorithm.

        Parameters
        ----------
        expression_data : pd.DataFrame
            Gene expression data stored in a data frame with sample identifiers as indices and 
            gene symbols as column names.

        Returns
        -------
        inferred_network : pd.DataFrame
            A data frame with gene symbols as indices and column names whose entries correspond to
            edge scores in inferred network.
        """

        # remove columns with zero standard deviation and normalize columns to unit variance
        expression_data = expression_data.loc[:, (expression_data.std() != 0)]
        expression_data = prp.normalizeToUnitVariance(expression_data)

        # get regulators and remove from list of regulators such genes that are not present in expression_data
        gene_names = np.array(expression_data.columns.copy())
        regulators = pd.read_csv('http://humantfs.ccbr.utoronto.ca/download/v_1.01/TFs_Ensembl_v_1.01.txt', sep='\t')
        regulators = np.array(regulators)[np.isin(regulators, gene_names)]

        targets = gene_names.tolist()
        regulators = regulators.tolist()
        if len(regulators) < 1 or len(gene_names) < 1:
            print('no overlap of regulators and genes in data set. No results for current block.')
    
        # save expression_data as csv
        main = os.path.join(test_suite, '..')
        data_path = os.path.join(main, 'temp', 'genie3_expression_data.csv') # TODO make unpredictable
        expression_data = expression_data.T # R version of genie wants gene x sample data set
        expression_data.to_csv(data_path, sep='\t')

        # GENIE3 read.expr.matrix TODO: set parameter nTrees
        out_path = os.path.join(main, 'temp', 'genie3_link_list.csv')
        genie3 = 'Rscript ../algorithms/GENIE3_R/GENIE3_script.R' 
        command = f'{genie3} {data_path} {out_path} {regulators} {targets}'
        subprocess.call(command, shell=True, stdout=subprocess.PIPE) # TODO something is not working

        # get results
        network = pd.read_csv(out_path, sep='\t')

        # remove temporary files
        subprocess.call('rm '+str(output_path), shell=True)
        subprocess.call('rm '+str(data_path), shell=True)

        return network

    def _get_top_k_edges(self, i, k):
            """Abstract method to return the top k edges for the inferred network for block i. 
            Must be implemented by derived classes.
            
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
                For undirected edges, ensure that <gene_2> <= <gene_2> for all tuples contained in edge set.
            """

            if self._inferred_networks is None:
                print('Call method infer_networks first.')
                return
            elif len(self.partition) <= i:
                print('Block with index ' + str(i) + ' does not exist in partition.')
                return

            block = self._inferred_networks[i]
            if block.shape[1] > 2:
                block = block.iloc[:, :2]
            if k >= len(block):
                k = len(block) - 1

            return set(block.iloc[:k, :])
