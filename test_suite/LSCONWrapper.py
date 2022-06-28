from NetworkInferenceWrapper import NetworkInferenceWrapper
import pandas as pd
import numpy as np
import os

#https://bitbucket.org/sonnhammergrni/lscon/src/master/
#https://academic.oup.com/bioinformatics/article/38/8/2263/6530276
#https://www.researchgate.net/profile/Daniel_Morgan6/publication/316465473_GeneSPIDER_-_Gene_regulatory_network_inference_benchmarking_with_controlled_network_and_data_properties/links/59c390830f7e9b21a82fcbf2/GeneSPIDER-Gene-regulatory-network-inference-benchmarking-with-controlled-network-and-data-properties.pdf?_sg%5B0%5D=tsz6yyqIVAxUrJa_wozjQ5uLeL2oc-ENf1AQEXjyNSNHna_iHNnlTFi_GbzdK-vT1tdnjYl3On3YURhY4FDTRA.VTTHtRSKiznFAb-wwl_lmUE3x4Fw1EVsd8jhDvyGXxezRq_znLWyEkSYBW50ELO-TGwSAYy1wL0TjWrGdgK4cA&_sg%5B1%5D=SPMlq3A9WnU7Qs3TOuQpna6Wj5kQuYm7GeAfzK4fCq5YK7lHCQA9fqecagXtvjg8oMnsv1xFmlupuBgU97fmm_u1H0cOXqnljci8U-dEoo2G.VTTHtRSKiznFAb-wwl_lmUE3x4Fw1EVsd8jhDvyGXxezRq_znLWyEkSYBW50ELO-TGwSAYy1wL0TjWrGdgK4cA&_iepl=
class LSCONWrapper(NetworkInferenceWrapper):

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
        prefix = 'lscon'
        out_path = os.path.join(main, 'temp', f'{prefix}_link_list.csv')

        # data

        # net

        # zetavec?

        # rawZeta?

        # regPath

        # run lscon
        cur = os.getcwd()
        os.chdir(os.path.join(main, 'algorithms', 'lscon'))
        command = None
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
            For undirected edges, ensure that <gene_1> <= <gene_2> for all tuples contained in edge set.
        """
        pass
