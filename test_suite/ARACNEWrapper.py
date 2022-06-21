from NetworkInferenceWrapper import NetworkInferenceWrapper
import pandas as pd
import numpy as np

#https://github.com/califano-lab/ARACNe-AP
#https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-7-S1-S7#Sec13
class ARACNEWrapper(NetworkInferenceWrapper):

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


        # -e: prepare .txt file with genes on rows and samples on columns; tab separated, with row- and column names
        # TODO: check again if normalization is done by the tool
        expression_data = expression_data.T
        expression_path = None # TODO

        # -o: set output folder
        main = os.path.join(test_suite, '..')
        prefix = 'aracne'
        out_path = os.path.join(main, 'temp', f'{prefix}_link_list.csv')

        # --tfs: prepare .txt file with known transcription factors and hand path
        tfs = pd.read_csv('http://humantfs.ccbr.utoronto.ca/download/v_1.01/TFs_Ensembl_v_1.01.txt', sep='\t')
        tfs_path = None # TODO

        # --pvalue: 1E-8 (from tutorial)
        p = 1E-8

        # run aracne: java -Xmx5G -jar aracne.jar -e test/matrix.txt  -o outputFolder --tfs test/tfs.txt --pvalue 1E-8
        cur = os.getcwd()
        os.chdir(os.path.join(main, 'algorithms', 'ARACNe-AP'))
        # TODO: copied from amim
        command = f'java -Xmx2g -jar keypathwayminer-standalone-5.0.jar -e {expression_path} -tfs {tfs_path} -o {out_path} --pvalue {p}'
        subprocess.run(command, shell=True)
        os.chdir(cur)

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
            k = min(k, len(block))
            top_k_edges = []
            for j in range(k):
                top_k_edges.append((block.iloc[j, 0], block.iloc[j, 1]))

            return set(top_k_edges)
