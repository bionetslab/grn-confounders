from abc import ABC, abstractmethod
import itertools as itt
import os

class NetworkInferenceWrapper(ABC):
    """An abstract wrapper class for network inference methods.
    Attributes
    ----------
    expression_data : pd.DataFrame
        Gene expression data stored in a data frame with sample identifiers as indices and 
        gene symbols as column names. Must be set explicitly before calling infer_networks().
    partition : list
        A list of lists containing sample identifiers for samples that fall into the different 
        blocks of the sample partition. Must be set explicitly before calling infer_networks().
    _inferred_networks : dict
        A dict containing data frames with gene symbols as indices and column names whose entries correspond to
        edge scores in inferred network for the different blocks contained in self.partition. Keys are block identifiers.
    Methods
    -------
    infer_networks():
        Infers all networks for the stored partition and expression data.
    mean_jaccard_index_at_k(k):
        Returns the mean Jaccard index for the top k edges in the inferred networks.
    _get_top_k_edges(i, k):
        Abstract method to returns the top k edges for the inferred network for block i. 
        Must be implemented by derived classes.
    Static methods
    --------------
    _infer_network(expression_data):
        Abstract method to infer a network from expression data. Must be implemented by derived classes.
    """
    
    def __init__(self, cwd):
        """Constructor."""
        self.expression_data = None
        self.partition = None
        self._inferred_networks = dict()
        self.cwd = cwd
        
    def infer_networks(self, rank=0):
        """Infers all networks for the stored partition and expression data. Before, delete any network stored in self._inferred_networks.
        Parameters
        ----------
        rank : int
            Rank of executing process. Default is 0.
        """
        self._inferred_networks.clear()
        for block_id in self.partition.keys():
            block = self.partition[block_id]
            self._inferred_networks.update({block_id : self._infer_network(self.expression_data.loc[block], rank)})

    def mean_jaccard_index_at_k(self, k):
        """Returns the mean Jaccard index for the top k edges in the inferred networks.
        Parameters
        ----------
        k : int
            Number of edges which should be considered when computing Jaccard indices.
        Return
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
        print(self._inferred_networks.keys())
        for i, j in itt.combinations(self._inferred_networks.keys(), 2):
            top_k_edges_i, state_i = self._get_top_k_edges(i, k)
            top_k_edges_j, state_j = self._get_top_k_edges(j, k)
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
    
    @staticmethod
    @abstractmethod        
    def _infer_network(expression_data, rank):
        """Abstract method to infer a network from expression data. Must be implemented by derived classes.
        Parameters
        ----------
        expression_data : pd.DataFrame
            Gene expression data stored in a data frame with sample identifiers as indices and 
            gene symbols as column names.
        Returns
        -------
        inferred_network : pd.DataFrame
            A data frame whose entries in the column 'score' correspond to edge scores in inferred network. An edge connects the genes
            whose symbols are in the columns 'source' and 'target', respectively, of the same row. The last column, 'type', indicates
            whether the network is directed or undirected.
        """
        pass
        
    @abstractmethod
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
            For undirected edges, ensure that <gene_1> <= <gene_2> for all tuples contained in edge set.
        """
        pass


            


