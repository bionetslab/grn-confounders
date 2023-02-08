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

    def __init__(self, data_dict, conf_dict, params_dict, rank=0):
        """Constructs TestRunner object on the requested parameters.
        
        Parameters
        ----------
        cwd : str
            Current working directory containing all input and output directories.
        data_dict : dict
            Dict of cancer types to be investigated with corresponding gene expression file names, phenotype file names, pheno data 
            filtering information, and separator character.
            e.g. {'tcga': False, 'ged': 'ged_filename.csv', 'pt': 'pt_filename.csv', 'sep': ',', 'tissue_type_field': 'sample_type.sample', 'tissue_type': 'Primary Tumor'}. Further documentation on the data_dict can be found at TODO.
            Important: Place ged and pt in data/ in root directory. 'tissue_type_field' must be a field in pt. Only samples with 'tissue_type' in 'tissue_type_field' will remain in the data set, i.e. will be used in the tests.
        conf_dict : dict
            Dict of confounders to be investigated with corresponding block types.
            e.g. {'gender.demographic': {'role': 'confounder', 'type': 'CATEGORY'}} Further documentation on the conf_dict can be found at TODO.
        params_dict : dict
            Dict of further parameters needed to run the tests.
            e.g. {'algorithms': [GENIE3, WGCNA], 'N_from': 0, 'N_to': 100, 'M_from': 0, 'M_to': 10, 'k_max': 5000, 'combine': False, 'par': False, 'g_all': False, 'save_networks': False} Further documentation on the params_dict can be found at TODO.
        """
        # set utility parameters
        self.cwd = os.getcwd()
        self.rank = rank
        self.save = params_dict['save_networks']
        self.n_from = params_dict['N_from']
        self.n_to = params_dict['N_to']
        self.m_from = params_dict['M_from']
        self.m_to = params_dict['M_to']
        self.k_max = params_dict['k_max']
        self.g_all = params_dict['g_all']
        self.combine = params_dict['combine']
        self.data_dict = data_dict
        self.conf_dict = conf_dict

        # required k in step size of 100, starting at 10
        self.rep_k = range(10, self.k_max, 100)

        # if g_all flag is set, add 'ALL' to confounders
        #if self.g_all: # TODO remove this. ALL is still needed, if so wants to infer entire network, but this line is not needed anymore
            #self.conf_dict[Selectors.ConfounderSelector.ALL] = {'role': Selectors.Role.CONFOUNDER, 'type': Selectors.BlockType.ALL}
        
        self.chi = [key for key in list(conf_dict.keys()) if conf_dict[key]['role'] == Selectors.Role.VARIABLE]
 
        # initialize and empty logfile
        self.logger = TestRunner.get_logger(params_dict['logfile'], self.cwd)

        # set selector parameters
        self.cancer_type_selectors = sorted(list(data_dict.keys()))
        self.confounder_selectors = list(conf_dict.keys())
        self.algorithm_selectors = [Selectors.AlgorithmSelector(val) for val in params_dict['algorithms']]

        # log start of test runs
        self._log_init()

        # prepare container for g_all networks
        self.g_all_networks = {ct_sel: {alg_sel: {} for alg_sel in self.algorithm_selectors} for ct_sel in self.cancer_type_selectors}

        # get data
        self.expression_datasets = {sel: self.get_expression_data(sel, data_dict[sel], self.logger) for sel in self.cancer_type_selectors}
        self.pheno_datasets = {sel: self.get_pheno_data(sel, data_dict[sel], self.logger) for sel in self.cancer_type_selectors}
        self.algorithm_wrappers = {sel: Selectors.get_algorithm_wrapper(sel) for sel in self.algorithm_selectors}

        # align data, remove 0-std genes, add combined data, if required
        self.expression_datasets, self.pheno_datasets = self.align_data(self.expression_datasets, self.pheno_datasets)
        self._add_if_combined()

        # prepare partition containers
        self.conf_partitions = {ct_sel: {conf_sel: '' for conf_sel in self.confounder_selectors} for ct_sel in self.cancer_type_selectors}
        self.rnd_partitions = {ct_sel: {conf_sel: '' for conf_sel in self.confounder_selectors} for ct_sel in self.cancer_type_selectors}
        
        # initialize container dicts for results
        self.conf_results = {ct_sel: {conf_sel: {alg_sel: {j: list([]) for j in range(self.m_from, self.m_to)} 
            for alg_sel in self.algorithm_selectors} for conf_sel in self.confounder_selectors} for ct_sel in self.cancer_type_selectors}
        self.rnd_results = {ct_sel: {conf_sel: {alg_sel: {i: list([]) for i in range(self.n_from, self.n_to)} 
            for alg_sel in self.algorithm_selectors} for conf_sel in self.confounder_selectors} for ct_sel in self.cancer_type_selectors}

    def add_custom_algorithm(self, wrapper, name=Selectors.AlgorithmSelector.CUSTOMWRAPPER):
        self.algorithm_selectors.append(name)
        self.algorithm_wrappers.update({name: wrapper})
        [alg.update({name: {j: list([]) for j in range(self.m_from, self.m_to)}}) for ct in self.conf_results.values() 
                    for conf in ct.values() for alg in conf.values()]
        [alg.update({name: {i: list([]) for i in range(self.n_from, self.n_to)}}) for ct in self.rnd_results.values() 
                    for conf in ct.values() for alg in conf.values()]

    def induce_partitions(self):
        """Induce partitions for all tests specified in the TestRunner object. If G_ALL is set, G_all is only run on confounder partition, 
        since confounder partition = random partition of the entire data."""
        self.conf_partitions = {ct_sel: {conf_sel: OrderedDict({ret[0]: ret[1] for ret in self.get_conf_partition(self.pheno_datasets[ct_sel], self.conf_dict[conf_sel]['type'], conf_sel, self.rank, logger=self.logger)}) 
            for conf_sel in self.confounder_selectors} for ct_sel in self.cancer_type_selectors}
        self.rnd_partitions = {ct_sel: {conf_sel: [OrderedDict({ret[0]: ret[1] for ret in self.get_ith_random_partition(i, self.pheno_datasets[ct_sel], self.conf_partitions[ct_sel][conf_sel], ct_sel, conf_sel)}) for i in range(self.n_from, self.n_to)]
            for conf_sel in self.confounder_selectors} for ct_sel in self.cancer_type_selectors}

    def infer_g_all(self):
        """If self.g_all is True, infer G_all for all cohorts specified in self.cancer_type_selectors with all algorithms
        specified in self.algorithm_selectors. Infer repeatedly for h in range(max(self.n_from, self.m_from), \
        min(self.n_to, self.m_to)) and save network in self.g_all_networks[ct_sel][alg_sel][h].
        """
        if not os.path.exists(os.path.join(self.cwd, 'results')):
            os.mkdir(os.path.join(self.cwd, 'results'))
        if not os.path.exists(os.path.join(self.cwd, 'results', 'JI')):
            os.mkdir(os.path.join(self.cwd, 'results', 'JI'))
        for ct_sel in self.cancer_type_selectors:
            for alg_sel in self.algorithm_selectors:
                if max(self.n_from, self.m_from) < min(self.n_to, self.m_to):
                    for alg_sel in self.algorithm_selectors:
                        print(f'Inferring G_all for cohort = {str(ct_sel)} with algorithm = {str(alg_sel)}')
                        algorithm_wrapper = self.algorithm_wrappers[alg_sel]
                        algorithm_wrapper.expression_data = self.expression_datasets[ct_sel]
                        algorithm_wrapper.partition = OrderedDict({ret[0]: ret[1] for ret in self.get_conf_partition(self.pheno_datasets[ct_sel], Selectors.BlockType.ALL, Selectors.ConfounderSelector.ALL, self.rank, self.logger)})
                        for h in range(max(self.n_from, self.m_from), min(self.n_to, self.m_to)):
                            algorithm_wrapper.infer_networks(self.rank)
                            self.save_networks(algorithm_wrapper._inferred_networks, h, 'conf', alg_sel, ct_sel, Selectors.ConfounderSelector.ALL, self.save)
                            self.g_all_networks[ct_sel][alg_sel].update({h: algorithm_wrapper._inferred_networks['all']})
                else:
                    self.logger.info('Comparison with G_all cannot be made for non-overlapping partition indices. Specify from, to such \
                        that max(self.n_from, self.m_from) < min(self.n_to, self.m_to).')

    def run_all(self):
        """Runs the tests for a given cancer_type and confounder on all algorithms. Test is only performed if the partition induced
        by a confounder contains more tahn one block. Performs all requested chi2 tests.
        """
        if not os.path.exists(os.path.join(self.cwd, 'results')):
            os.mkdir(os.path.join(self.cwd, 'results'))
        if not os.path.exists(os.path.join(self.cwd, 'results', 'JI')):
            os.mkdir(os.path.join(self.cwd, 'results', 'JI'))
        if self.g_all:
            self.infer_g_all()
        for ct_sel in self.cancer_type_selectors:
            for conf_sel in self.confounder_selectors:
                if len(list(self.conf_partitions[ct_sel][conf_sel].values())) > 1 or conf_sel == Selectors.ConfounderSelector.ALL:
                    self._run_chi2_tests(conf_sel, ct_sel)
                    self._run_on_cancer_type_confounder(ct_sel, conf_sel)

    def _run_on_cancer_type_confounder(self, ct_sel, conf_sel):
        """Runs the tests for a given cancer type and confounder with all specified algorithms.
        Parameters
        ----------
        ct_sel : str | CancerTypeSelector
            Cohort that schoulb be investigated.
        conf_sel : str | ConfounderSelector
            Confounder that should be investigated.
        """
        for alg_sel in self.algorithm_selectors:
            algorithm_wrapper = self.algorithm_wrappers[alg_sel]
            algorithm_wrapper.expression_data = self.expression_datasets[ct_sel]
            print(f'\t\talgorithm = {str(alg_sel)}')
            print('running on confounder-based partitions...')
            algorithm_wrapper.partition = self.conf_partitions[ct_sel][conf_sel]
            for j in range(self.m_from, self.m_to):
                algorithm_wrapper.infer_networks(self.rank)
                self.save_networks(algorithm_wrapper._inferred_networks, j, 'conf', alg_sel, ct_sel, conf_sel, self.save)
                network_state = []
                intersections = []
                unions = []
                if conf_sel == Selectors.ConfounderSelector.ALL: # TODO remove? 
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
                
                if self.g_all and j in self.g_all_networks[ct_sel][alg_sel].keys():
                    cur_g_all = self.g_all_networks[ct_sel][alg_sel][j]
                    block_networks = algorithm_wrapper._inferred_networks.copy()
                    for block_id in block_networks.keys():
                        algorithm_wrapper._inferred_networks = {'g_all': cur_g_all, block_id: block_networks[block_id]}
                        # compute meanJIs over k of cur_block_network and cur_g_all
                        network_state = []
                        intersections = []
                        unions = []
                        results = []
                        for k in self.rep_k:
                            ji, state, s_int, s_un = algorithm_wrapper.mean_jaccard_index_at_k(k)
                            results.append(ji)
                            unions.append(s_un)
                            intersections.append(s_int)
                            network_state.append(state)
                        pd.DataFrame({'size intersection': intersections, 'size union': unions, 'state': network_state, 'k': self.rep_k, 
                        'mean JI': results}).to_csv(os.path.join('results', 'JI', f'g_all_conf_{str(j)}_{str(alg_sel)}_{str(conf_sel)}_{str(ct_sel)}_{block_id}_jaccInd.csv'), index=False)
                        
            print('running on random partitions...')
            if conf_sel != Selectors.ConfounderSelector.ALL:
                for i in range(self.n_from, self.n_to):
                    algorithm_wrapper.partition = self.rnd_partitions[ct_sel][conf_sel][i-self.n_from]
                    algorithm_wrapper.infer_networks(self.rank)
                    self.save_networks(algorithm_wrapper._inferred_networks, i, 'rnd', alg_sel, ct_sel, conf_sel, self.save)
                    network_state = []
                    intersections = []
                    unions = []
                    for k in self.rep_k:
                        ji, state, s_int, s_un = algorithm_wrapper.mean_jaccard_index_at_k(k)
                        self.rnd_results[ct_sel][conf_sel][alg_sel][i].append(ji)
                        unions.append(s_un)
                        intersections.append(s_int)
                        network_state.append(state)
                    pd.DataFrame({'size intersection': intersections, 'size union': unions, 'state': network_state, 'k': self.rep_k, 
                    'mean JI': self.rnd_results[ct_sel][conf_sel][alg_sel][i]}).to_csv(os.path.join('results', 'JI', 
                    f'rnd_{i}_{str(alg_sel)}_{str(conf_sel)}_{str(ct_sel)}_jaccInd.csv'), index=False)
    
                    if self.g_all and i in self.g_all_networks[ct_sel][alg_sel].keys():
                        cur_g_all = self.g_all_networks[ct_sel][alg_sel][j]
                        block_networks = algorithm_wrapper._inferred_networks.copy()
                        for block_id in block_networks.keys():
                            algorithm_wrapper._inferred_networks = {'g_all': cur_g_all, block_id: block_networks[block_id]}
                            # compute meanJIs over k of cur_block_network and cur_g_all
                            network_state = []
                            intersections = []
                            unions = []
                            results = []
                            for k in self.rep_k:
                                ji, state, s_int, s_un = algorithm_wrapper.mean_jaccard_index_at_k(k)
                                results.append(ji)
                                unions.append(s_un)
                                intersections.append(s_int)
                                network_state.append(state)
                            pd.DataFrame({'size intersection': intersections, 'size union': unions, 'state': network_state, 'k': self.rep_k, 'mean JI': results}).to_csv(os.path.join('results', 'JI', f'g_all_rnd_{str(i)}_{str(alg_sel)}_{str(conf_sel)}_{str(ct_sel)}_{block_id}_jaccInd.csv'), index=False)
                            
    def save_networks(self, inferred_networks, part_nb, mode, alg_sel, ct_sel, conf_sel, save=False):
        """Saves the inferred networks to csv.
        Parameters
        ----------
        inferred_networks: dict
            dict containing the inferred networks as pd.DataFrames. Key is the identifier of the block that was used for network inference.
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
        save : bool
            Whether to save or not save the inferred networks. The networks are usually very large and saving all networks might not be necessary and/or possible. Defaults to False.
        """
        if not os.path.exists(os.path.join(self.cwd, 'results')):
            os.mkdir(os.path.join(self.cwd, 'results'))
        if not os.path.exists(os.path.join(self.cwd, 'results', 'networks')):
            os.mkdir(os.path.join(self.cwd, 'results', 'networks'))
        if save:
            for block_id in inferred_networks.keys():
                path = os.path.join(self.cwd, 'results', 'networks', f'{mode}_part{part_nb}_{block_id}_{alg_sel}_{ct_sel}_{conf_sel}_gene_list.csv')
                top_k_edges = inferred_networks[block_id].iloc[:self.k_max, :]
                top_k_edges.to_csv(path, index = False)

    def get_expression_data(self, cancer_type_selector, sel_dict, logger=None):
        """Loads the expression data for the selected cancer type. Only leaves columns of protein-coding genes, provided
        in the file protein-coding_gene.csv, in expression_data. Removes gene version identifiers from gene ensembl IDs.
        Parameters
        ----------
        cancer_type_selector : str
            Specifies for which cancer type the phenotypes should be loaded.
        sel_dict : dict
            Dictionary containing information about at least ged (expression data file name) and sep (separator used in ged file).
        logger : logging.Logger
            Logger to be used to log data specs and progress.
        Returns
        -------
        expression_data : pd.DataFrame
            Expression data (indices are sample IDs, column names are gene IDs).
        """
        if not os.path.exists(os.path.join(self.cwd, 'partitions')):
            os.mkdir(os.path.join(self.cwd, 'partitions'))
        expression_data = pd.read_csv(os.path.join(self.cwd, 'data', sel_dict['ged']), sep=sel_dict['sep'], header=0, index_col=0)
        print('Remove version identifiers from gene symbols in expression data for cohort ' + str(cancer_type_selector) + '...')
        expression_data.columns = expression_data.columns.str.split('.').str[0].tolist()
        print('Only leave protein-coding genes in expression data set for cohort ' + str(cancer_type_selector) + '...')
        pcg = pd.read_csv(os.path.join('data', 'protein-coding_gene.csv'))
        mask = expression_data.columns.intersection(pcg['ensembl_gene_id'].values)
        expression_data = expression_data[mask]
        if logger:
            logger.info(cancer_type_selector + ' - expression data from file: ' + sel_dict['ged'] + ' - data shape after removing non-coding genes: ')
            logger.info(expression_data.shape)
            logger.info('\n')
        return expression_data

    def get_pheno_data(self, cancer_type_selector, sel_dict, logger=None):
        """Loads the phenotype data for the selected cancer type and adds a column containing cancer_type_selector.
        Parameters
        ----------
        cancer_type_selector : str
            Specifies for which cancer type the phenotypes should be loaded.
        sel_dict : dict
            Dictionary containing information about at least pt (pheno data file name) and sep (separator used in ged file).
            Optional tissue_type_field and tissue_type to filter pheno data.
        logger : logging.Logger
            Logger to be used to log data specs and progress.
        Returns
        -------
        pheno_data : pd.DataFrame
            Expression data (indices are sample IDs, column names are gene IDs).
        """
        try:
            tissue_type_field, tissue_type = sel_dict['tissue_type_field'], sel_dict['tissue_type']
        except:
            tissue_type_field, tissue_type = None, None
        if logger:
            logger.info(cancer_type_selector + ' - Reading pheno type data from file: ' + sel_dict['pt'])
        pheno_data = pd.read_csv(os.path.join(self.cwd, 'data', sel_dict['pt']), sep=sel_dict['sep'], header=0, index_col=0)
        assert len(pheno_data.iloc[0]) == len(pheno_data.iloc[0].values)
        pheno_data['cohort'] = str(cancer_type_selector)
        if tissue_type is not None and tissue_type_field is not None:
            try:
                if logger:
                    logger.info('Filter pheno type data for ' + tissue_type + ' accrding to field ' + tissue_type_field + '\n')
                pheno_data =  pheno_data[pheno_data[tissue_type_field] == tissue_type]
            except KeyError:
                print('Filtering failed: ' + tissue_type_field + ' not in pheno data. Continue with unfiltered pheno data.')
                if logger:
                    logger.info('Filtering failed: ' + tissue_type_field + ' not in pheno data. Continue with unfiltered data.\n')
        if logger:
            logger.info('Pheno data shape: ')
            logger.info(str(pheno_data.shape) + '\n')
        return pheno_data

    def get_conf_partition(self, pheno_data_orig, block_type, pheno_field, rank=0, min_block_size=20, logger=None):
        """Returns two lists with the first containing string-identifiers for the blocks of the requested confounder 
        and the second containing the sample ids corresponding to the blocks. For the age confounder, the lower and upper quartiles are
        computed separately per cohort and are then combined into the resulting upper and lower fragments.
        Parameters
        ----------
        pheno_data_orig : pd.DataFrame
            Data frame containing phenotypic information. One row per sample, one column per attribute.
        block_type : Selectors.BlockType
            QUARTILE or CATEGORY; defines how to create the partition.
        pheno_field : str
            Field in pheno type file to be used to induce the partition.
        rank : int
            Rank of the executing process. Default is 0.
        min_block_size : int
            Minimum block size. If a block is smaller than min_block_size, it is removed from the partition.
        logger : logging.Logger
            Logger to be used to log data specs and progress.
        Returns
        -------
        conf_partition : list
            Contains the blocks belonging to the confounder-based partition.
        """
        if logger:
            logger.info('Induce partition by ' + str(pheno_field) +'\n')
        pheno_data = pheno_data_orig.copy()
        indices = None
        blocks = []
        conf_partition = []
        if block_type == Selectors.BlockType.ALL or pheno_field == str(Selectors.ConfounderSelector.ALL):
            samples = pheno_data.index.tolist()
            conf_partition.append(('all', samples))
            if logger:
                logger.info('Do not create blocks, but use entire data\n')
            return conf_partition
        pheno_data = pheno_data[pheno_data[pheno_field] != 'not reported']
        pheno_data = pheno_data[pheno_data[pheno_field].notna()]
        if block_type == Selectors.BlockType.CATEGORY:
            blocks = sorted(list(set(pheno_data[pheno_field].str.strip().values)))
            if logger:
                logger.info('Induce partition by ' + str(pheno_field) + '\n')
            for block_attr in blocks:
                samples = pheno_data.loc[pheno_data[pheno_field].str.strip() == block_attr].index.tolist()
                if len(samples) >= min_block_size:
                    conf_partition.append((block_attr, samples))
                    if logger:
                        logger.info('block ' + block_attr + ': ' + str(len(samples)) + ' samples\n')
        elif block_type == Selectors.BlockType.QUARTILE:
            samples_lower = []
            samples_upper = []
            for cohort in set(pheno_data['cohort'].str.strip().values):
                pheno_cohort = pheno_data[pheno_data['cohort'] == cohort]
                lower, upper = pheno_cohort[pheno_field].quantile(0.25), pheno_cohort[pheno_field].quantile(0.75)
                samples_lower.extend(pheno_cohort.loc[pheno_cohort[pheno_field] <= lower].index.tolist())
                samples_upper.extend(pheno_cohort.loc[pheno_cohort[pheno_field] > upper].index.tolist())
            if len(samples_lower) >= min_block_size and len(samples_upper) >= min_block_size:
                conf_partition.append(('lower', samples_lower))
                conf_partition.append(('upper', samples_upper))
                if logger:
                    logger.info('block lower: ' + str(len(samples_lower)) + ' samples\n')
                    logger.info('block upper: ' + str(len(samples_upper)) + ' samples\n')
        return conf_partition

    def get_ith_random_partition(self, i, samples, conf_partition_dict, ct_sel, conf_sel):
        """Returns n random partitions each containing blocks of the same size as in the corresponding
        confounder based partition.
        Parameters
        ----------
        i : int
            Specifies the index of the random partition to be generated and/or returned.
        samples : pd.DataFrame
            Pheno data set of the cohort to be partitioned.
        conf_partition_dict : dict
            List of blocks as pd.DataFrames with one column containing the sample identifiers belonging to the block.
        ct_sel : str
            String identifier of cancer type (cohort).
        conf_sel : str
            String identifier of confounder.
        Returns
        -------
        partition : list
            List of tuples with block identifiers and randomly sampled block. Block sizes according to blocks of the confounder-based partition.
        """
        samples_cpy = samples.copy()
        cur = []
        try:
            if not os.path.exists(os.path.join(self.cwd, 'partitions')):
                os.mkdir(os.path.join(self.cwd, 'partitions'))
            part = pd.read_csv(os.path.join(self.cwd, 'partitions', f'rnd_part{i}_{ct_sel}_{conf_sel}'), header=None, index_col=False, dtype=str).values.tolist()
            begin = 0
            end = 0
            conf_partition = list(conf_partition_dict.values())
            for key in conf_partition_dict.keys():
                end += len(conf_partition_dict[key])
                block = [item for sublist in part[begin:end] for item in sublist]
                cur.append((key, block))
                begin += len(conf_partition_dict[key])
        except FileNotFoundError:
            print(f'rnd_partition {i} not found. Create new partitions.')
            for key in conf_partition_dict.keys():
                block = samples_cpy.sample(n=len(conf_partition_dict[key]), replace=False).index.values
                samples_cpy = samples_cpy.drop(block)
                cur.append((key, block))
                pd.DataFrame(block).to_csv(os.path.join(self.cwd, 'partitions', f'rnd_part{i}_{ct_sel}_{conf_sel}'), mode='a', header=False, index=False)
        
        return cur

    def align_data(self, expression_datasets, pheno_datasets, logger=None):
        """Data alignment. Remove such samples from the expression_data files that are not in the pheno_data files and vice versa.
        Remove all genes where standard deviation is 0.
        
        Parameters
        ----------
        expression_datasets : pd.DataFrame
            Expression data set to be aligned by samples with pheno_datasets.

        pheno_datasets : pd.DataFrame
            Pheno data set to be aligned by samples with expression_datasets.

        Returns
        -------
        expression_datasets : pd.DataFrame
            Expression data set aligned by samples with pheno_datasets.

        pheno_datasets : pd.DataFrame
            Pheno data set aligned by samples with expression_datasets.
        """
        for sel in expression_datasets.keys():
            print('Align expression data and phenotype data on samples for cohort ' + str(sel) + '...')
            pheno_datasets[sel] = pheno_datasets[sel][pheno_datasets[sel].index.isin(expression_datasets[sel].index)]
            samples = pheno_datasets[sel].index.values        
            expression_datasets[sel] = expression_datasets[sel].loc[samples]
            print('Remove genes where standard deviation of expression data is 0 for cohort ' + str(sel) + '...')
            expression_datasets[sel] = expression_datasets[sel].loc[:, (expression_datasets[sel].std() != 0)]
            pheno_datasets[sel] = pheno_datasets[sel][pheno_datasets[sel].index.isin(expression_datasets[sel].index)]
            if logger:
                logger.info('Pheno data shape of ' + str(sel) + ' after alignment: ' + str(pheno_datasets[sel].shape) + '\nGene expression data shape of ' + str(sel) + ' after alignment: ' + str(expression_datasets[sel].shape))
        return expression_datasets, pheno_datasets

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
                p, table = self._var_conf_chi(self.pheno_datasets[ct_sel], conf_sel, var, self.conf_dict, self.logger)

    def _var_conf_chi(self, pheno, conf_sel, var, conf_dict, logger=None):
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
        p = 10000000
        if logger:
            logger.info('Chi^2 test of ' + str(conf_sel) + ' and ' + str(var) + '\n')
        conf_partition = self.get_conf_partition(pheno, conf_dict[conf_sel]['type'], conf_sel)
        confusion_table = pd.DataFrame()
        for block in conf_partition:
            var_partition = self.get_conf_partition(pheno.loc[block[1]], conf_dict[var]['type'], var, min_block_size=0)
            for var_block in var_partition:
                confusion_table.loc[var_block[0], block[0]] = len(var_block[1])
        try:
            confusion_table = confusion_table.dropna()
            chi2, p, dof, ex = chi2_contingency(confusion_table, correction=False)
            sign = ('**' if p < 0.01 else '*') if p < 0.05 else 'ns'
            if logger:
                logger.info(str(confusion_table) + '\n')
                logger.info(['chi^2 test ' + str(conf_sel) +' and ' + var + ' significant: ' + str(p) + f'{sign}\n'])
        except:
            if logger:
                logger.info('No chi^2 test of ' + str(conf_sel) + ' and ' + str(var) + ' possible.\n')
        return p, confusion_table

    @staticmethod
    def get_logger(logfile, cwd):
        """Initialize and return logger that prints at INFO level to a file in cwd named logfile.
        Parameters
        ----------
        logfile : str
            Name of the file to log to. Is placed into cwd.
        cwd : str
            Current working directory containing all input and output directories.
        Return
        ----------
        logger : logging.Logger
            INFO logger that prints to logfile.
        """
        logger = logging.getLogger(logfile)
        handler = logging.FileHandler(os.path.join(cwd, logfile), mode='w')
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
            self.logger.info(str(ct) + '\n')
        self.logger.info('Confounders:\n')
        for conf in self.confounder_selectors:
            self.logger.info(str(conf) + '\n')
        self.logger.info('Algorithms:\n')
        for alg in self.algorithm_selectors:
            self.logger.info(str(alg) + '\n')
