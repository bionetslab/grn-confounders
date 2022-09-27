from . import Selectors
import pandas as pd
import numpy as np
import os

class TestRunner(object):
    """Runs the tests."""

    def __init__(self, cancer_types, confounders, algorithms, n_from, n_to, m_from, m_to, k, combine, rank):
        """Constructs TestRunner object on the requested parameters.
        
        Parameters
        ----------
        cancer_types : list
            List of strings that specify the cohorts that should be investigated.

        confounders : list
            List of strings that specify the confounders that should be investigated.

        n_from : int
            Starting index for the TestRunner.

        n_to : int
            Ending index for the TestRunner.

        k : int
            Number of edges to compare in the resulting networks.

        combine : bool
            Indicates whether to run all tests on the combination of all datasets, plus an additional test on the cancer type
            variable.

        rank : int
            Rank identifier of the process.
        """
        self.rank = rank
        self.n_from = n_from
        self.n_to = n_to
        self.m_from = m_from
        self.m_to = m_to
        self.k = k

        self.cancer_type_selectors = [Selectors.CancerTypeSelector(val) for val in cancer_types]
        self.confounder_selectors = [Selectors.ConfounderSelector(val) for val in confounders]
        self.algorithm_selectors = [Selectors.AlgorithmSelector(val) for val in algorithms]

        self.expression_datasets = {sel: Selectors.get_expression_data(sel) for sel in self.cancer_type_selectors}
        self.pheno_datasets = {sel: Selectors.get_pheno_data(sel) for sel in self.cancer_type_selectors}
        self.algorithm_wrappers = {sel: Selectors.get_algorithm_wrapper(sel) for sel in self.algorithm_selectors}

        if combine:
            comb_name = '-'.join([str(el) for el in self.cancer_type_selectors])
            self.confounder_selectors.append(Selectors.ConfounderSelector.TYPE)
            self.cancer_type_selectors.append(comb_name)
            self.expression_datasets.update({comb_name: pd.concat(self.expression_datasets.values())})
            self.pheno_datasets.update({comb_name: pd.concat(self.pheno_datasets.values())})

        self.preprocessData()
        
        self.conf_partitions = {ct_sel: {conf_sel: Selectors.get_conf_partition(self.pheno_datasets[ct_sel], conf_sel) for conf_sel in self.confounder_selectors}
            for ct_sel in self.cancer_type_selectors}
        self.rnd_partitions = {ct_sel: {conf_sel: Selectors.get_n_random_partitions(self.n_from, self.n_to, self.pheno_datasets[ct_sel]['submitter_id.samples'], self.conf_partitions[ct_sel][conf_sel], ct_sel, conf_sel)
            for conf_sel in self.confounder_selectors} for ct_sel in self.cancer_type_selectors}
        
        self.conf_results = {ct_sel: {conf_sel: {alg_sel: {j: list([]) for j in range(self.m_from, self.m_to)} for alg_sel in self.algorithm_selectors} for conf_sel in self.confounder_selectors} for ct_sel in self.cancer_type_selectors} # a list of length n per cancer_type and confounder
        self.rnd_results = {ct_sel: {conf_sel: {alg_sel: {i: list([]) for i in range(self.n_from, self.n_to)} for alg_sel in self.algorithm_selectors} for conf_sel in self.confounder_selectors} for ct_sel in self.cancer_type_selectors}
        
    def run_on_cancer_types_confounders(self, combine, verbose):
        """Runs the tests for a given cancer_type and confounder on all algorithms. Test is only performed if the partition induced
        by a confounder contains more tahn one block.

        Parameters
        ----------
        combine : bool
            Perform tests on combination of all cancer type datasets and additionally, run the tests for the variable "cancer type" as a confounder.

        verbose : bool
            Print progress to stdout.
        """
        for ct_sel in self.cancer_type_selectors:
            for conf_sel in self.confounder_selectors:
                if len(self.rnd_partitions[ct_sel][conf_sel][0]) > 1:
                    self.run_on_all_cancer_types_confounders_partitions(ct_sel, conf_sel, combine, verbose)

    def run_on_all_cancer_types_confounders_partitions(self, ct_sel, conf_sel, combine, verbose=False):
        """Runs the tests for a given cancer type and confounder with all specified algorithms.

        Parameters
        ----------
        ct_sel : CancerTypeSelector
            Cohort that schoulb be investigated.

        conf_sel : ConfounderSelector
            Confounder that should be investigated.

        combine : bool
            Perform tests on combination of all cancer type datasets and additionally, run the tests for the variable "cancer type" as a confounder.

        verbose : bool
            Print progress to stdout.
        """
        for alg_sel in self.algorithm_selectors:
            prefix = f'{str(alg_sel)}'
            algorithm_wrapper = self.algorithm_wrappers[alg_sel]
            algorithm_wrapper.expression_data = self.expression_datasets[ct_sel]
            if verbose:
                print(f'\t\talgorithm = {str(alg_sel)}')

            if verbose:
                print('running on random partitions...')
            for i in range(self.n_from, self.n_to):
                algorithm_wrapper.partition = self.rnd_partitions[ct_sel][conf_sel][i-self.n_from]
                algorithm_wrapper.infer_networks(self.rank)
                self.save_networks(algorithm_wrapper._inferred_networks, i, 'rnd', alg_sel, ct_sel, conf_sel)
                index = []
                network_state = []
                intersections = []
                unions = []
                for k in range(10, self.k, 10):
                    ji, state, s_int, s_un = algorithm_wrapper.mean_jaccard_index_at_k(k)
                    self.rnd_results[ct_sel][conf_sel][alg_sel][i].append(ji)
                    index.append(k)
                    unions.append(s_un)
                    intersections.append(s_int)
                    network_state.append(state)
                pd.DataFrame({'size intersection': intersections, 'size union': unions, 'state': network_state, 'k': index, 'mean JI': self.rnd_results[ct_sel][conf_sel][alg_sel][i]}).to_csv(os.path.join('results', 'JI', f'rnd_{i}_{str(alg_sel)}_{str(conf_sel)}_{str(ct_sel)}_jaccInd.csv'), index=False)
            if verbose:
                print('running on confounder-based partitions...')
            algorithm_wrapper.partition = self.conf_partitions[ct_sel][conf_sel]
            for j in range(self.m_from, self.m_to):
                algorithm_wrapper.infer_networks(self.rank)
                self.save_networks(algorithm_wrapper._inferred_networks, j, 'conf', alg_sel, ct_sel, conf_sel)
                index=[]
                network_state = []
                intersections = []
                unions = []
                for k in range(10, self.k, 50):
                    ji, state, s_int, s_un = algorithm_wrapper.mean_jaccard_index_at_k(k)
                    self.conf_results[ct_sel][conf_sel][alg_sel][j].append(ji)
                    index.append(k)
                    unions.append(s_un)
                    intersections.append(s_int)
                    network_state.append(state)
                pd.DataFrame({'size intersection': intersections, 'size union': unions, 'state': network_state, 'k': index, 'mean JI': self.conf_results[ct_sel][conf_sel][alg_sel][j]}).to_csv(os.path.join('results', 'JI', f'cb_{j}_{str(alg_sel)}_{str(conf_sel)}_{str(ct_sel)}_jaccInd.csv'), index=False)
                
            if combine:
                if verbose:
                    print('starting tests on the combination of all datasets...')
                self.run_cancer_type_as_confounder()

    def run_cancer_type_as_confounder(self):
        pass

    def preprocessData(self):
        """Data preprocessing. Remove such samples from the expression_data files that are not in the pheno_data files and vice versa. Removes all 
        samples of type other than 'Primary Tumor from pheno_data. Removes version identifiers from the gene symbols in expression_data."""
        for sel in self.cancer_type_selectors:
            print('Align expression data and phenotype data on samples for cohort ' + str(sel) + '...')
            keep = self.pheno_datasets[sel]['submitter_id.samples'].isin(self.expression_datasets[sel].index)
            self.pheno_datasets[sel] = self.pheno_datasets[sel][keep]
            self.pheno_datasets[sel] = self.pheno_datasets[sel][self.pheno_datasets[sel]['submitter_id.samples'].isin(self.expression_datasets[sel].index)]
            samples = self.pheno_datasets[sel]['submitter_id.samples']            
            self.expression_datasets[sel] = self.expression_datasets[sel].loc[samples]
            print('Remove genes where standard deviation of expression data is 0 for cohort ' + str(sel) + '...')
            self.expression_datasets[sel] = self.expression_datasets[sel].loc[:, (self.expression_datasets[sel].std() != 0)]
            self.pheno_datasets[sel] = self.pheno_datasets[sel][self.pheno_datasets[sel]['submitter_id.samples'].isin(self.expression_datasets[sel].index)]

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
        """Clears the results of the previous run."""
        self.rnd_results = None
        self.conf_results = None
