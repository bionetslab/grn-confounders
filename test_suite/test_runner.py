from . import Selectors
import pandas as pd
import numpy as np
import os

class TestRunner(object):
    """Runs the tests."""

    def __init__(self, n_from, n_to, m_from, m_to, k):
        """Constructs TestRunner object."""
        self.n_from = n_from
        self.n_to = n_to
        self.m_from = m_from
        self.m_to = m_to
        self.k = k # max k to select in top_k_edges
        
        self.cancer_type_selectors = list(Selectors.CancerTypeSelector) # a "cohort" in TCGA terminology # TODO: rename to 'cohort'?
        self.algorithm_selectors = list(Selectors.AlgorithmSelector)
        self.confounder_selectors = list(Selectors.ConfounderSelector)
        
        self.expression_datasets = {sel: Selectors.get_expression_data(sel) for sel in self.cancer_type_selectors}
        self.pheno_datasets = {sel: Selectors.get_pheno_data(sel) for sel in self.cancer_type_selectors}
        self.preprocessData()
        self.algorithm_wrappers = {sel: Selectors.get_algorithm_wrapper(sel) for sel in self.algorithm_selectors}
        self.conf_partitions = {ct_sel: {conf_sel: Selectors.get_conf_partition(self.pheno_datasets[ct_sel], conf_sel) for conf_sel in self.confounder_selectors}
            for ct_sel in self.cancer_type_selectors}
        self.rnd_partitions = {ct_sel: {conf_sel: Selectors.get_n_random_partitions(self.n_from, self.n_to, self.pheno_datasets[ct_sel]['submitter_id.samples'], self.conf_partitions[ct_sel][conf_sel], ct_sel, conf_sel)
            for conf_sel in self.confounder_selectors} for ct_sel in self.cancer_type_selectors}
        
        self.cancer_type_names = []
        self.algorithm_names = []
        self.confounder_names = []
        self.conf_results = {ct_sel: {conf_sel: {alg_sel: list([]) for alg_sel in self.algorithm_selectors} for conf_sel in self.confounder_selectors} for ct_sel in self.cancer_type_selectors} # a list of length n per cancer_type and confounder
        self.rnd_results = {ct_sel: {conf_sel: {alg_sel: {i: list([]) for i in range(self.n_from, self.n_to)} for alg_sel in self.algorithm_selectors} for conf_sel in self.confounder_selectors} for ct_sel in self.cancer_type_selectors}
        self.outfile = ''
        
    def run_all(self):
        """Runs the tests for all cancer types and all confounders on all algorithms.
        """
        for ct_sel in self.cancer_type_selectors:
            self.cancer_type_names.append(str(ct_sel))
            for conf_sel in self.confounder_selectors:
                self.confounder_names.append(str(conf_sel))
                self.run_on_all_cancer_types_confounders_partitions(ct_sel, conf_sel, verbose)

    def run_on_cancer_types_confounders(self, cancer_types, confounders, algorithms, verbose):
        """Runs the tests for a given cancer_type and confounder on all algorithms.

        Parameters
        ----------
        cancer_types : list
            List of strings that specify the cohorts that should be investigated.

        confounders : list
            List of strings that specify the confounders that should be investigated.

        algorithms: list
            List of strings that specify the algorithms that should be investigated.

        verbose : bool
            Print progress to stdout.
        """
        ct_selectors = [Selectors.CancerTypeSelector(val) for val in cancer_types]
        conf_selectors = [Selectors.ConfounderSelector(val) for val in confounders]
        for ct_sel in ct_selectors:
            self.cancer_type_names.append(str(ct_sel))
            for conf_sel in conf_selectors:
                self.confounder_names.append(str(conf_sel))
                self.run_on_all_cancer_types_confounders_partitions(ct_sel, conf_sel, algorithms, verbose)

    def run_on_all_cancer_types_confounders_partitions(self, ct_sel, conf_sel, algorithms, verbose=False):
        """Runs the tests for a given cancer_type and confounder on all algorithm.

        Parameters
        ----------
        ct_sel : CancerTypeSelector
            Cohort that schoulb be investigated.

        conf_sel : ConfounderSelector
            Confounder that should be investigated.
            
        algorithms: list
            List of strings that specify the algorithms that should be investigated.

        verbose : bool
            Print progress to stdout.
        """
        alg_selectors = [Selectors.AlgorithmSelector(val) for val in algorithms]
        for alg_sel in alg_selectors:
            self.algorithm_names.append(str(alg_sel))
            prefix = f'{str(alg_sel)}'
            if verbose:
                print(f'\t\talgorithm = {str(alg_sel)}')
            
            algorithm_wrapper = self.algorithm_wrappers[alg_sel]
            #if str(alg_sel) == 'GENIE3':
                #algorithm_wrapper.expression_data = self.expression_datasets[ct_sel].iloc[:, :5000]
            #else:
            algorithm_wrapper.expression_data = self.expression_datasets[ct_sel]
            
            print('running on random partitions...')
            for i in range(self.n_from, self.n_to):
                algorithm_wrapper.partition = self.rnd_partitions[ct_sel][conf_sel][i]
                algorithm_wrapper.infer_networks()
                self.save_networks(algorithm_wrapper._inferred_networks, i, 'rnd', alg_sel, ct_sel, conf_sel)
                index = []
                for k in range(10, self.k, 10):
                    try:
                        self.rnd_results[ct_sel][conf_sel][alg_sel][i].append(algorithm_wrapper.mean_jaccard_index_at_k(k))
                        index.append(k)
                    except IndexError:
                        self.rnd_results[ct_sel][conf_sel][alg_sel][i].append(-1)
                        index.append(k)
                pd.DataFrame({'k': index, 'mean JI': self.rnd_results[ct_sel][conf_sel][alg_sel][i]}).to_csv(os.path.join('results', 'JI', f'rnd_{i}_{str(alg_sel)}_{str(conf_sel)}_{str(ct_sel)}_jaccInd.csv'), index=False)
        
            print('running on confounder-based partitions...')
            algorithm_wrapper.partition = self.conf_partitions[ct_sel][conf_sel]
            for j in range(self.m_from, self.m_to):
                algorithm_wrapper.infer_networks()
                self.save_networks(algorithm_wrapper._inferred_networks, 0, 'conf', alg_sel, ct_sel, conf_sel)
                index=[]
                for k in range(10, self.k, 10):
                    try:
                        self.conf_results[ct_sel][conf_sel][alg_sel].append(algorithm_wrapper.mean_jaccard_index_at_k(k))
                        index.append(k)
                    except IndexError:
                        self.rnd_results[ct_sel][conf_sel][alg_sel][i].append(-1)
                        index.append(k)
                pd.DataFrame({'k': index, 'mean JI': self.conf_results[ct_sel][conf_sel][alg_sel]}).to_csv(os.path.join('results', 'JI', f'cb_{j}_{str(alg_sel)}_{str(conf_sel)}_{str(ct_sel)}_jaccInd.csv'), index=False)
                
    def preprocessData(self):
        """Data preprocessing. Remove such samples from the expression_data files that are not in the pheno_data files and vice versa. Removes all 
        samples of type other than 'Primary Tumor from pheno_data. Removes version identifiers from the gene symbols in expression_data."""
        for sel in self.cancer_type_selectors:
            print('Filter Primary Tumor samples in pheno data for cohort ' + str(sel) + '...')
            self.pheno_datasets[sel] =  self.pheno_datasets[sel][self.pheno_datasets[sel]['sample_type.samples'] == 'Primary Tumor']
            print('Remove version identifiers from gene symbols in expression data for cohort ' + str(sel) + '...')
            self.expression_datasets[sel].columns = self.expression_datasets[sel].columns.str.split('.').str[0].tolist()
            print('Align expression data and phenotype data on samples for cohort ' + str(sel) + '...')
            self.expression_datasets[sel] = self.expression_datasets[sel].loc[self.expression_datasets[sel].index.intersection(self.pheno_datasets[sel]['submitter_id.samples'])]
            self.pheno_datasets[sel] = self.pheno_datasets[sel][self.pheno_datasets[sel]['submitter_id.samples'].isin(self.expression_datasets[sel].index)]
            self.expression_datasets[sel] = self.expression_datasets[sel].loc[:, (self.expression_datasets[sel].std() != 0)]

    def save_networks(self, inferred_networks, part_nb, mode, alg_sel, ct_sel, conf_sel):
        """Saves the inferred networks to csv.

        Parameters
        ----------
        inferred_networks: list
            list containing the inferred networks as pd.DataFrames

        part_nb: int
            index of the random partition whose results are saved; must be between 0 and self.n. If confounder-based partition, part_nb is 0.

        mode: str
            'rnd' for random partition, 'conf' for confounder-based partition.

        alg_sel: AlgorithmSelector
            algorithm that produced the results

        ct_sel : CancerTypeSelector
            Cohort that was investigated.

        conf_sel : ConfounderSelector
            Confounder that was investigated.
        """
        print('saving the results')
        cwd = os.getcwd()
        for block_nb in range(len(inferred_networks)):
            path = os.path.join(cwd, 'results', 'networks', f'{mode}_part{part_nb}_block{block_nb}_{alg_sel}_{ct_sel}_{conf_sel}_gene_list.csv')
            inferred_networks[block_nb].to_csv(path, index = False)

    def clear(self):
        """Clears the results of the last previous run."""
        self.cancer_type_names = []
        self.algorithm_names = []
        self.confounder_names = []
        self.rnd_results = None
        self.conf_results = None
        self.outfile = ''
