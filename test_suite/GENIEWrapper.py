from NetworkInferenceWrapper import NetworkInferenceWrapper

# NOTE: we want to provide a script that works specifically on TCGA data; it takes the data and puts it into the format that
# is required here. As this framework should also be applicable to data from other sources, this preprocessing procedure
# is not included here. This preprocessing script should also perform the removal of version ids from the gene ids.
# TODO: ask DB about his opinion

# TODO: define path name conventions

class GENIEWrapper(NetworkInferenceWrapper):

    def _infer_network(expression_data):

        """Abstract method to infer a network from expression data. Must be implemented by derived classes.
        
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


        # 1. bring gene expression data into right format, get gene- and sample identifiers and create pd.DataFrame df
            # note: from now on, we only cut the data frame to make sure that the order of the sample names and gene names 
            # always match the data
            # remove 0-std cols


        # 2. get known-TFs as regulators


        # 3. bring phenotype data into right format


        # 4. remove non-"Primary Tumor" samples from df


        # 5. form confounder-based partitions
            # normalize to unit varince
            # remove 0-std columns


        # 6. form random partitions


        # 7. call GENIE3 on both confounder-based partitions and random partitions


        # 8. save results in permanent file


    def _get_top_k_edges(self, i, k):
            """Abstract method to returns the top k edges for the inferred network for block i. 
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