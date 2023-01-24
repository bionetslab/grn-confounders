from . import Selectors
from . import preprocessing
import pandas as pd
import numpy as np
import os
from scipy.stats import chi2_contingency
from datetime import datetime
from collections import OrderedDict
import logging

class TestRunner(object):
    """Runs the tests."""

    def __init__(self, data_dict, conf_dict, algorithms, n_from, n_to, m_from, m_to, k, tissue_type_field, tissue_type, combine=False, g_all=False, rank=0, logfile='logfile.txt'):
        """Constructs TestRunner object on the requested parameters.
        
        Parameters
        ----------
        data_dict : dict
            Dict of cancer types to be investigated with corresponding gene expression file names, phenotype file names, and separator character.
        conf_dict : dict
            Dict of confounders to be investigated with corresponding block types.
        n_from : int
            Starting index for the TestRunner.
        n_to : int
            Ending index for the TestRunner.
        k : int
            Number of edges to compare in the resulting networks.
        sep : str
            Character to be used as separator for files. Default is \',\'.
        tissue_type_field : str
            Field in pt file to be used to filter the data. Default is None.
        tissue_type : str
            Attribute to be filtered for in the tissue_type_field. Default is None.
        logfile : str
            Name of file to use as logfile. Opened in mode w at first log.
        combine : bool
            Indicates whether to run all tests on the combination of all datasets, plus an additional test on the cancer type
            variable.
        rank : int
            Rank identifier of the process.
        g_all : bool
            Boolean indicating whether a network should be inferred from the entire gene expression data or not.
        """
        # set utility parameters
        self.rank = rank
        self.n_from = n_from
        self.n_to = n_to
        self.m_from = m_from
        self.m_to = m_to
        self.k = k
        self.g_all = g_all
        self.combine = combine
        self.data_dict = data_dict
        self.conf_dict = conf_dict

        # required k in step size of 100, starting at 10
        self.rep_k = range(10, self.k, 100)

        # if g_all flag is set, add 'NONE' to confounders
        if self.g_all:
            self.conf_dict[Selectors.ConfounderSelector.NONE] = Selectors.BlockType.ALL
        self.chi = [key for key in conf_dict.keys() if conf_dict[key]['role'] == 'variable']

        # initialize and empty logfile
        self.logger = TestRunner.get_logger(logfile)

        # set selector parameters
        self.cancer_type_selectors = sorted(list(data_dict.keys()))
        self.confounder_selectors = list(conf_dict.keys())
        self.algorithm_selectors = [Selectors.AlgorithmSelector(val) for val in algorithms]

        self._log_init()

        # get data
        self.expression_datasets = {sel: Selectors.get_expression_data(sel, data_dict[sel]['ged'], data_dict[sel]['sep'], self.logger) for sel in self.cancer_type_selectors}
        self.pheno_datasets = {sel: Selectors.get_pheno_data(sel, data_dict[sel]['pt'], data_dict[sel]['sep'], tissue_type_field, tissue_type, self.logger) for sel in self.cancer_type_selectors}
        self.algorithm_wrappers = {sel: Selectors.get_algorithm_wrapper(sel) for sel in self.algorithm_selectors}

        # align data, remove 0-std genes, add combined data, if required
        self.expression_datasets, self.pheno_datasets = Selectors.align_data(self.expression_datasets, self.pheno_datasets)
        self._add_if_combined()

        # induce partitions for all tests. If G_ALL is set, G_all is only run on confounder partition, since confounder partition
        # = random partition entire data
        self.conf_partitions = {ct_sel: {conf_sel: OrderedDict({ret[0]: ret[1] for ret in Selectors.get_conf_partition(self.pheno_datasets[ct_sel], conf_dict[conf_sel]['type'], conf_sel, self.rank, self.logger)}) 
            for conf_sel in self.confounder_selectors} for ct_sel in self.cancer_type_selectors}
        self.rnd_partitions = {ct_sel: {conf_sel: Selectors.get_n_random_partitions(self.n_from, self.n_to, self.pheno_datasets[ct_sel], list(self.conf_partitions[ct_sel][conf_sel].values()), ct_sel, conf_sel)
            for conf_sel in self.confounder_selectors} for ct_sel in self.cancer_type_selectors}
        
        # initialize container dicts for results
        self.conf_results = {ct_sel: {conf_sel: {alg_sel: {j: list([]) for j in range(self.m_from, self.m_to)} 
            for alg_sel in self.algorithm_selectors} for conf_sel in self.confounder_selectors} for ct_sel in self.cancer_type_selectors}
        self.rnd_results = {ct_sel: {conf_sel: {alg_sel: {i: list([]) for i in range(self.n_from, self.n_to)} 
            for alg_sel in self.algorithm_selectors} for conf_sel in self.confounder_selectors} for ct_sel in self.cancer_type_selectors}

    def run_on_cancer_types_confounders(self):
        """Runs the tests for a given cancer_type and confounder on all algorithms. Test is only performed if the partition induced
        by a confounder contains more tahn one block. Performs all requested chi2 tests.
        """
        for ct_sel in self.cancer_type_selectors:
            for conf_sel in self.confounder_selectors:
                if len(list(self.conf_partitions[ct_sel][conf_sel].values())) > 1 or conf_sel == Selectors.ConfounderSelector.NONE:
                    self._run_chi2_tests(conf_sel, ct_sel)
                    self.run_on_all_cancer_types_confounders_partitions(ct_sel, conf_sel)

    def run_on_all_cancer_types_confounders_partitions(self, ct_sel, conf_sel):
        """Runs the tests for a given cancer type and confounder with all specified algorithms.
        Parameters
        ----------
        ct_sel : str | CancerTypeSelector
            Cohort that schoulb be investigated.
        conf_sel : str | ConfounderSelector
            Confounder that should be investigated.
        """
        for alg_sel in self.algorithm_selectors:
            prefix = f'{str(alg_sel)}'

            algorithm_wrapper = self.algorithm_wrappers[alg_sel]
            algorithm_wrapper.expression_data = self.expression_datasets[ct_sel]
            print(f'\t\talgorithm = {str(alg_sel)}')

            print('running on confounder-based partitions...')
            algorithm_wrapper.partition = list(self.conf_partitions[ct_sel][conf_sel].values())
            for j in range(self.m_from, self.m_to):
                algorithm_wrapper.infer_networks(self.rank)
                self.save_networks(algorithm_wrapper._inferred_networks, j, 'conf', alg_sel, ct_sel, conf_sel)
                network_state = []
                intersections = []
                unions = []
                if conf_sel == Selectors.ConfounderSelector.NONE:
                    continue
                for k in self.rep_k:
                    ji, state, s_int, s_un = algorithm_wrapper.mean_jaccard_index_at_k(k)
                    self.conf_results[ct_sel][conf_sel][alg_sel][j].append(ji)
                    unions.append(s_un)
                    intersections.append(s_int)
                    network_state.append(state)
                pd.DataFrame({'size intersection': intersections, 'size union': unions, 'state': network_state, 'k': self.rep_k, 
                'mean JI': self.conf_results[ct_sel][conf_sel][alg_sel][j]}).to_csv(os.path.join('results', 'JI', 
                f'cb_{j}_{str(alg_sel)}_{str(conf_sel)}_{str(ct_sel)}_jaccInd.csv'), index=False)

            print('running on random partitions...')
            if conf_sel != Selectors.ConfounderSelector.NONE:
                for i in range(self.n_from, self.n_to):
                    algorithm_wrapper.partition = self.rnd_partitions[ct_sel][conf_sel][i-self.n_from]
                    algorithm_wrapper.infer_networks(self.rank)
                    self.save_networks(algorithm_wrapper._inferred_networks, i, 'rnd', alg_sel, ct_sel, conf_sel)
                    network_state = []
                    intersections = []
                    unions = []
                    # skip random partition for G_all, since for G_all, random partition = confounder partition
                    for k in self.rep_k:
                        ji, state, s_int, s_un = algorithm_wrapper.mean_jaccard_index_at_k(k)
                        self.rnd_results[ct_sel][conf_sel][alg_sel][i].append(ji)
                        unions.append(s_un)
                        intersections.append(s_int)
                        network_state.append(state)
                    pd.DataFrame({'size intersection': intersections, 'size union': unions, 'state': network_state, 'k': self.rep_k, 
                    'mean JI': self.rnd_results[ct_sel][conf_sel][alg_sel][i]}).to_csv(os.path.join('results', 'JI', 
                    f'rnd_{i}_{str(alg_sel)}_{str(conf_sel)}_{str(ct_sel)}_jaccInd.csv'), index=False)
    
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
        ct_sel : str | CancerTypeSelector
            Cohort that was investigated.
        conf_sel : str | ConfounderSelector
            Confounder that was investigated.
        """
        cwd = os.getcwd()
        for block_nb in range(len(inferred_networks)):
            block_id = str(list(self.conf_partitions[ct_sel][conf_sel].keys())[block_nb])
            path = os.path.join(cwd, 'results', 'networks', f'{mode}_part{part_nb}_{block_id}_{alg_sel}_{ct_sel}_{conf_sel}_gene_list.csv')
            inferred_networks[block_nb].to_csv(path, index = False)

    def _add_if_combined(self):
        """Combine cohorts: concat gene expression data and phenotype data and add TYPE confounder to test effect of cohort type."""
        if self.combine:
            comb_name = '-'.join([str(el) for el in self.cancer_type_selectors])
            self.conf_dict[comb_name] = Selectors.BlockType.CATEGORY
            self.confounder_selectors.append(Selectors.ConfounderSelector.TYPE)
            self.cancer_type_selectors.append(comb_name)
            self.expression_datasets.update({comb_name: pd.concat(self.expression_datasets.values())})
            self.pheno_datasets.update({comb_name: pd.concat(self.pheno_datasets.values())})

    def _run_chi2_tests(self, conf_sel, ct_sel):
        if self.rank == 0:
            for var in self.chi:
                p, table = TestRunner._var_conf_chi(self.pheno_datasets[ct_sel], conf_sel, var, self.conf_dict)

    @staticmethod
    def _var_conf_chi(pheno, conf_sel, var, conf_dict):
        """Perform Chi^2 test of conf_sel and var.
        Parameters
        ----------
        pheno : pd.DataFrame
            Phenotype file of the cohort to test.
        conf_sel : str
            confounder to test. Must correspond to a field (column name) in pheno.
        var : str
            variable to test. Must correspond to a field (column name) in pheno.
        conf_dict : dict
            Dictionary with at least two keys var and conf_sel, and a nested dict with at least key \'type\' containing a
            Selectors.BlockType.
        Return
        ----------
        p : float
            P-value of the conducted Chi^2 test.
        confusion_table : pd.DataFrame
            Initial confusion table of the conducted Chi^2 test.
        """
        if self.logger:
            self.logger.info('Chi^2 test of ' + str(conf_sel) + ' and ' + str(var) + '\n')
        conf_partition = Selectors.get_conf_partition(pheno, conf_dict[conf_sel]['type'], conf_sel)
        confusion_table = pd.DataFrame()
        for block in conf_partition:
            var_partition = Selectors.get_conf_partition(pheno.loc[block[1]], conf_dict[var]['type'], var, min_block_size=0)
            for var_block in var_partition:
                confusion_table.loc[var_block[0], block[0]] = len(var_block[1])
        try:
            confusion_table = confusion_table.dropna()
            chi2, p, dof, ex = chi2_contingency(confusion_table, correction=False)
            sign = ('**' if p < 0.01 else '*') if p < 0.05 else 'ns'
            if self.logger:
                self.logger.info(str(confusion_table) + '\n')
                self.logger.info(['chi^2 test ' + str(conf_sel) +' and ' + var + ' significant: ' + str(p) + f'{sign}\n'])
        except:
            if self.logger:
                self.logger.info('No chi^2 test of ' + str(conf_sel) + ' and ' + str(var) + ' possible.\n')
        return p, confusion_table

    @staticmethod
    def get_logger(logfile):
        """Initialize and return logger that prints at INFO level to a file in cwd named logfile.
        Parameters
        ----------
        logfile : str
            Name of the file to log to. Is placed into cwd.
        Return
        ----------
        logger : logging.Logger
            INFO logger that prints to logfile.
        """
        logger = logging.getLogger(logfile)
        handler = logging.FileHandler(os.path.join(os.getcwd(), logfile), mode='w')
        handler.setLevel(logging.INFO)
        fmt = logging.Formatter('%(message)s')
        handler.setFormatter(fmt)
        logger.addHandler(handler)
        logger.setLevel(logging.INFO)
        return logger

    def _log_init(self):
        """Log start date, cohorts, confounders, and algorithm selectors of current test run to logfile."""
        if not self.logger:
            print('No logger initialized. Will skip all logging.')
            return
        self.logger.info('Initialized test runner object - ')
        self.logger.info(datetime.now().strftime("%d/%m/%Y %H:%M:%S") + '\n')
        self.logger.info('Cohorts:\n')
        for ct in self.cancer_type_selectors:
            self.logger.info(ct + '\n')
        self.logger.info('Confounders:\n')
        for conf in self.confounder_selectors:
            self.logger.info(conf + '\n')
        self.logger.info('Algorithms:\n')
        for alg in self.algorithm_selectors:
            self.logger.info(alg + '\n')
