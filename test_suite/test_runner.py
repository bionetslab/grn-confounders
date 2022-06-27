from . import Selectors
import pandas as pd
import numpy as np

class TestRunner(object):
    """Runs the tests."""

    def __init__(self, n):
        """Constructs TestRunner object."""
        self.n = n
        self.cancer_type_selectors = list(Selectors.CancerTypeSelector) # a "cohort" in TCGA terminology
        self.algorithm_selectors = list(Selectors.AlgorithmSelector)
        self.confounder_selectors = list(Selectors.ConfounderSelector)

        self.expression_datasets = {sel: Selectors.get_expression_data(sel) for sel in self.cancer_type_selectors}
        self.pheno_datasets = {sel: Selectors.get_pheno_data(sel) for sel in self.cancer_type_selectors} # only hand over Primary Tumor samples
        self.algorithm_wrappers = {sel: Selectors.get_algorithm_wrapper(sel) for sel in self.algorithm_selectors}
        self.conf_partitions = {ct_sel: {conf_sel: Selectors.get_conf_partition(self.pheno_datasets[ct_sel], conf_sel) for conf_sel in self.confounder_selectors}
            for ct_sel in self.cancer_type_selectors}
        self.rnd_partitions = {ct_sel: [Selectors.get_n_random_partitions(self.n, self.pheno_datasets[ct_sel]['submitter_id.samples'], self.conf_partitions[ct_sel][conf_sel])
            for conf_sel in self.confounder_selectors] for ct_sel in self.cancer_type_selectors}

        self.cancer_type_names = []
        self.algorithm_names = []
        self.confounder_names = []
        self.partitions = []
        self.conf_partition = []
        self.results = {ct_sel: {conf_sel: [] for conf_sel in self.confounder_selectors} for ct_sel in self.cancer_type_selectors} # a list of length n per cancer_type and confounder
        self.outfile = ''

    def run_on_all_cancer_types_confounders_partitions(self, n, cancer_type_selector, confounder_selector, algorithm_selector, verbose):
        """Runs the tests for a given cancer_type, confounder and algorithm on n random partitions and the
        confounder based partition.

        Parameters
        ----------
        n : int
            Number of random partitions to be generated.
        cancer_type_selector : CancerTypeSelector
            Specifies the cancer type that should be investigated.
        confounder_selector : ConfounderSelector
            Specifies the confounder that should be investigated.
        algorithm_selector : AlgorithmSelector
            Specifies the algorithm that should be run.
        verbose : bool
            Print progress to stdout.
        """
        Selectors.download_known_tfs()

        self.cancer_type_names.append(str(cancer_type_selector))
        self.algorithm_names.append(str(algorithm_selector))
        self.confounder_names.append(str(confounder_selector))

        #self.conf_partition = self.conf_partitions[cancer_type_selector][confounder_selector]
        #samples = self.pheno_datasets[cancer_type_selector][0, :]
        #self.partitions = Selectors.get_n_random_partitions(n, samples, self.conf_partition)

        prefix = f'{str(algorithm_selector)}'
        if verbose:
            print(f'\t\talgorithm = {str(algorithm_selector)}')

        algorithm_wrapper = self.algorithm_wrappers[algorithm_selector]
        algorithm_wrapper.expression_data = self.expression_datasets[cancer_type_selector]

        for partition in self.partitions + self.conf_partition: # TODO how to distinguish between random and confounder based condition
            algorithm_wrapper.partition = partition
            algorithm_wrapper.infer_networks()
            self.results[cancer_type_selector][confounder_selector].append(algorithm_wrapper.mean_jaccard_index_at_k(5)) # TODO exemplarily only for 5

    def clear(self):
        """Clears the results of the last previous run."""
        self.cancer_type_names = []
        self.algorithm_names = []
        self.confounder_names = []
        self.partitions = []
        self.results = None
        self.outfile = ''
