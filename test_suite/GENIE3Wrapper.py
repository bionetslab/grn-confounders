from .NetworkInferenceWrapper import NetworkInferenceWrapper
import sys
import os
import pandas as pd
import numpy as np
import subprocess
from . import preprocessing as prp
import csv
test_suite = os.path.join(os.path.dirname(__file__))
sys.path.append(test_suite)

# TODO: define path name conventions
# TODO make temporary file names unpredictable (?)

class GENIE3Wrapper(NetworkInferenceWrapper):

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
        main = os.path.join(test_suite, '..')
        prefix = 'genie3'

        # remove columns with zero standard deviation and normalize columns to unit variance
        expression_data = expression_data.loc[:, (expression_data.std() != 0)]
        expression_data = prp.normalizeToUnitVariance(expression_data)

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

        # TODO: hand paths to R script (won't recognize the paths)
        os.chdir(os.path.join(main, 'algorithms', 'GENIE3_R'))
        command = f'Rscript GENIE3_R_wrapper.R {prefix}'# {regulator_path} {out_path}'
        ret = subprocess.run(command, shell=True)
        os.chdir(test_suite)

        # get results
        network = pd.read_csv(out_path, sep='\t')
        
        # remove temporary files
        subprocess.call('rm '+str(out_path), shell=True)
        subprocess.call('rm '+str(data_path), shell=True)
        subprocess.call('rm '+str(regulator_path), shell=True)
        
        return network # TODO ask about "gene symbols as indices"
        
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
            block = self._inferred_networks[i] # TODO change in case that gene symbols as indices
            k = min(k, len(block))
            top_k_edges = []
            for j in range(k):
                top_k_edges.append((block.iloc[j, 0], block.iloc[j, 1]))
            return set(top_k_edges)
