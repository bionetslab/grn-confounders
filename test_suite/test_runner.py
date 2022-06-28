from . import Selectors
import pandas as pd
import numpy as np

class TestRunner(object):
    """Runs the tests."""

    def __init__(self, n):
        """Constructs TestRunner object."""
        self.n = n
        self.cancer_type_selectors = list(Selectors.CancerTypeSelector) # a "cohort" in TCGA terminology # TODO: rename to 'cohort'?
        self.algorithm_selectors = list(Selectors.AlgorithmSelector)
        self.confounder_selectors = list(Selectors.ConfounderSelector)

        self.expression_datasets = {sel: Selectors.get_expression_data(sel) for sel in self.cancer_type_selectors} # version identifiers in gene symbols lready cut off
        self.pheno_datasets = {sel: Selectors.get_pheno_data(sel) for sel in self.cancer_type_selectors} # only hand over Primary Tumor samples
        self.algorithm_wrappers = {sel: Selectors.get_algorithm_wrapper(sel) for sel in self.algorithm_selectors}
        self.conf_partitions = {ct_sel: {conf_sel: Selectors.get_conf_partition(self.pheno_datasets[ct_sel], conf_sel) for conf_sel in self.confounder_selectors}
            for ct_sel in self.cancer_type_selectors}
        self.rnd_partitions = {ct_sel: Selectors.get_n_random_partitions(self.n, self.pheno_datasets[ct_sel]['submitter_id.samples'], self.conf_partitions[ct_sel][conf_sel])
            for conf_sel in self.confounder_selectors for ct_sel in self.cancer_type_selectors}

        self.cancer_type_names = []
        self.algorithm_names = []
        self.confounder_names = []
        self.conf_results = {ct_sel: {conf_sel: [] for conf_sel in self.confounder_selectors} for ct_sel in self.cancer_type_selectors} # a list of length n per cancer_type and confounder
        self.rnd_results = {ct_sel: {conf_sel: [] for conf_sel in self.confounder_selectors} for ct_sel in self.cancer_type_selectors}
        self.outfile = ''

    def run_on_all_cancer_types_confounders_partitions(self, verbose=False):
        """Runs the tests for a given cancer_type, confounder and algorithm on n random partitions and the
        confounder based partition.

        Parameters
        ----------
        n : int
            Number of random partitions to be generated.

        verbose : bool
            Print progress to stdout.
        """
        Selectors.download_known_tfs()

        for ct_sel in self.cancer_type_selectors:
            self.cancer_type_names.append(str(ct_sel))
            for conf_sel in self.confounder_selectors:
                self.confounder_names.append(str(conf_sel))
                for alg_sel in self.algorithm_selectors:
                    self.algorithm_names.append(str(alg_sel))
                    prefix = f'{str(alg_sel)}'
                    if verbose:
                        print(f'\t\talgorithm = {str(alg_sel)}')

                    algorithm_wrapper = self.algorithm_wrappers[alg_sel]
                    algorithm_wrapper.expression_data = self.expression_datasets[ct_sel].iloc[:150, :500] #TODO: remove later
                    for col in algorithm_wrapper.expression_data: # TODO: this should be done after data download later
                        algorithm_wrapper.expression_data.rename({col: col.split('.')[0]}, axis=1, inplace=True)

                    print('running on confounder-based partitions...') # TODO: what do we do, if one block is empty?
                    algorithm_wrapper.partition = self.conf_partitions[ct_sel][conf_sel] # TODO: what if the intersection of regulators and genes in expr_data is empty?
                    algorithm_wrapper.infer_networks()
                    self.conf_results[ct_sel][conf_sel].append(algorithm_wrapper.mean_jaccard_index_at_k(5)) # TODO exemplarily only for 5

                    print('running on random partitions...')
                    for ct_sel in self.cancer_type_selectors:
                        for i in range(self.n):
                            algorithm_wrapper.partition = self.rnd_partitions[ct_sel][i]
                            algorithm_wrapper.infer_networks()
                            self.rnd_results[ct_sel][conf_sel].append(algorithm_wrapper.mean_jaccard_index_at_k(5)) # TODO exemplarily only for 5
                    
                    print(self.conf_results[ct_sel][conf_sel])

    def clear(self):
        """Clears the results of the last previous run."""
        self.cancer_type_names = []
        self.algorithm_names = []
        self.confounder_names = []
        self.rnd_results = None
        self.conf_results = None
        self.outfile = ''

"""
        cancer_type_selector : CancerTypeSelector
            Specifies the cancer type that should be investigated.

        confounder_selector : ConfounderSelector
            Specifies the confounder that should be investigated.

        algorithm_selector : AlgorithmSelector
            Specifies the algorithm that should be run.
"""