import Selectors
import preprocessing
import pandas as pd
import numpy as np
import os
import itertools as itt
from scipy.stats import chi2_contingency
from scipy.stats import mannwhitneyu
from scipy.stats import wilcoxon
from datetime import datetime
from collections import OrderedDict

class TestRunner(object):
    """Runs the tests."""

    def __init__(self, data_dict, conf_dict, algorithms, n_from, n_to, m_from, m_to, k, combine=False, sep='\t', g_all=False, chi=['tumor_stage.diagnoses'], rank=0, logfile='logfile.txt'):
        """Constructs TestRunner object on the requested parameters.
        
        Parameters
        ----------
        data_dict : dict
            Dict of cancer types to be investigated with corresponding gene expression file names and phenotype file names.

        conf_dict : dict
            Dict of confounders to be investigated with corresponding block types.

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

        sep : str
            Character symbol specifying the separator used in the gene expression files and phenotype files.

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
        self.chi = chi

        # initialize and empty logfile
        self.logfile = logfile
        self.stamp = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        self._log([self.stamp], mode='w')

        # set selector parameters
        self.cancer_type_selectors = sorted(list(data_dict.keys()))
        self.confounder_selectors = list(conf_dict.keys())
        self.algorithm_selectors = [Selectors.AlgorithmSelector(val) for val in algorithms]

        # get data
        self.expression_datasets = {sel: Selectors.get_expression_data(sel, data_dict[sel]['ged'], sep) for sel in self.cancer_type_selectors}
        self.pheno_datasets = {sel: Selectors.get_pheno_data(sel, data_dict[sel]['pt'], sep) for sel in self.cancer_type_selectors}
        self.algorithm_wrappers = {sel: Selectors.get_algorithm_wrapper(sel) for sel in self.algorithm_selectors}

        # align data, remove 0-std genes, add combined data, if required
        self.expression_datasets, self.pheno_datasets = Selectors.align_data(self.expression_datasets, self.pheno_datasets)
        self._add_if_combined()

        # induce partitions for all tests. If G_ALL is set, G_all is only run on confounder partition, since confounder partition
        # = random partition = entire data
        self.conf_partitions = {ct_sel: {conf_sel: OrderedDict({ret[0]: ret[1] for ret in Selectors.get_conf_partition(self.pheno_datasets[ct_sel], conf_dict[conf_sel], conf_sel, self.rank)}) 
            for conf_sel in self.confounder_selectors} for ct_sel in self.cancer_type_selectors}
        self.rnd_partitions = {ct_sel: {conf_sel: Selectors.get_n_random_partitions(self.n_from, self.n_to, self.pheno_datasets[ct_sel].index, list(self.conf_partitions[ct_sel][conf_sel].values()), ct_sel, conf_sel)
            for conf_sel in self.confounder_selectors} for ct_sel in self.cancer_type_selectors}
        
        # initialize container dicts for results
        self.conf_results = {ct_sel: {conf_sel: {alg_sel: {j: list([]) for j in range(self.m_from, self.m_to)} 
            for alg_sel in self.algorithm_selectors} for conf_sel in self.confounder_selectors} for ct_sel in self.cancer_type_selectors}
        self.rnd_results = {ct_sel: {conf_sel: {alg_sel: {i: list([]) for i in range(self.n_from, self.n_to)} 
            for alg_sel in self.algorithm_selectors} for conf_sel in self.confounder_selectors} for ct_sel in self.cancer_type_selectors}
        
        # initilize data frame to store fraction of significant p-values in one-sided Mann-Whitney-U test
        self.mwu_fractions = pd.DataFrame(columns=['cohort', 'confounder/variable', 'method', 'frac'])

    def run_on_cancer_types_confounders(self):
        """Runs the tests for a given cancer_type and confounder on all algorithms. Test is only performed if the partition induced
        by a confounder contains more tahn one block.
        """
        for ct_sel in self.cancer_type_selectors:
            for conf_sel in self.confounder_selectors:
                if len(list(self.conf_partitions[ct_sel][conf_sel].values())) > 1 or conf_sel == Selectors.ConfounderSelector.NONE:
                    # chi^2 tests
                    self._run_chi2_tests(conf_sel, ct_sel)
                    self.run_on_all_cancer_types_confounders_partitions(ct_sel, conf_sel)

            for conf_sel in self.confounder_selectors:
                for alg_sl in self.algorithm_selectors:
                    if conf_sel != Selectors.ConfounderSelector.NONE:
                        print('Run one-sided Mann-Whitney-U test for each k...')
                        self.run_mwu_tests(self, ct_sel, conf_sel, alg_sel)
                    if conf_sel != Selectors.ConfounderSelector.NONE:
                        print('Run Wilcoxon tests for each k...')
                        self.run_wilcoxon_tests(ct_sel, conf_sel, alg_sel)
                        
            self.mwu_fractions.to_csv(os.path.join('results', self.stamp + '_mwu_fractions.csv'))
            self.visualize_mwus(ct_sel)
            if self.g_all:
                self.wilcox_fractions.to_csv(os.path.join('results', self.stamp + '_wilcox_fractions.csv'))

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
            for i in range(self.n_from, self.n_to):
                algorithm_wrapper.partition = self.rnd_partitions[ct_sel][conf_sel][i-self.n_from]
                algorithm_wrapper.infer_networks(self.rank)
                self.save_networks(algorithm_wrapper._inferred_networks, i, 'rnd', alg_sel, ct_sel, conf_sel)
                network_state = []
                intersections = []
                unions = []
                # skip random partition for G_all, since for G_all, random partition = confounder partition
                if conf_sel != Selectors.ConfounderSelector.NONE:
                    for k in self.rep_k:
                        ji, state, s_int, s_un = algorithm_wrapper.mean_jaccard_index_at_k(k)
                        self.rnd_results[ct_sel][conf_sel][alg_sel][i].append(ji)
                        unions.append(s_un)
                        intersections.append(s_int)
                        network_state.append(state)
                    pd.DataFrame({'size intersection': intersections, 'size union': unions, 'state': network_state, 'k': self.rep_k, 
                    'mean JI': self.rnd_results[ct_sel][conf_sel][alg_sel][i]}).to_csv(os.path.join('results', 'JI', 
                    f'rnd_{i}_{str(alg_sel)}_{str(conf_sel)}_{str(ct_sel)}_jaccInd.csv'), index=False)

    def run_mwu_tests(self, ct_sel, conf_sel, alg_sel):
        r = pd.concat([pd.DataFrame({'mean JI':rnd, 'k':self.rep_k}) for rnd in self.rnd_results[ct_sel][conf_sel][alg_sel].values()])
        c = pd.concat([pd.DataFrame({'mean JI':conf, 'cohort':ct_sel, 'confounder/variable':str(conf_sel), 'method':alg_sel, 'k':self.rep_k}) for conf in self.conf_results[ct_sel][conf_sel][alg_sel].values()])
        assert len(c) == len(r)
        significants = 0
        for k in self.rep_k:
            r_k = r[r['k'] == k]
            c_k = c[c['k'] == k]
            mwu = mannwhitneyu(c_k['mean JI'],r_k['mean JI'], alternative='less')
            if mwu.pvalue < 0.05:
                significants += 1
        self.mwu_fractions = pd.concat([self.mwu_fractions, pd.DataFrame({'cohort':ct_sel, 'confounder/variable':str(conf_sel), 
                            'method':str(alg_sel), 'frac':significants/len(self.rep_k)})])
    
    def run_wilcoxon_tests(self, ct_sel, conf_sel, alg_sel):
        block_ids = list(self.conf_partitions[ct_sel][conf_sel].keys())
        print('Compute mean JIs for G_all and confounder-based networks...')
        networks_blocks_conf = {bl: pd.DataFrame(columns=['size intersection', 'size union', 'state', 'k', 'mean JI', 'partID']) for bl in block_ids}
        conf_results = {l:{bl:[] for bl in block_ids} for l in range(self.m_from, self.m_to)}
        for j in range(self.m_from, self.m_to):
            plain = pd.read_csv(os.path.join(os.getcwd(), 'results', 'networks', 'conf_part'+str(j)+'_blockall_'+alg_sel+'_'+ct_sel+'_none_gene_list.csv'), header=0)
            for bl in list(self.conf_partitions[ct_sel][conf_sel].keys()):
                block_network = os.path.join(os.getcwd(), 'results', 'networks', 'conf_part'+str(j)+'_block'+bl+'_'+alg_sel+'_'+ct_sel+'_'+conf_sel+'_gene_list.csv')
                inferred_networks = [plain, pd.read_csv(block_network, header=0)]
                for k in self.rep_k:
                    ji, state, s_int, s_un = mean_jaccard_index_at_k(k, inferred_networks, alg_sel)
                    conf_results[j][bl].append(ji)
                JI_conf = pd.DataFrame({'k': self.rep_k, 'mean JI': conf_results[j][bl]}, 'partID':j,'confounder':conf_sel, 
                                        'method':alg_sel, 'cohort':ct_sel)
                networks_blocks_conf[bl] = pd.concat(networks_blocks_conf[bl], JI_conf)

        print('Compute mean JIs for G_all and random networks...')
        networks_blocks_rnd = {bl: pd.DataFrame(columns=['size intersection', 'size union', 'state', 'k', 'mean JI', 'partID']) for bl in block_ids}
        rnd_results = {l:{bl:[] for bl in block_ids} for l in range(self.n_from, self.n_to)}
        for j in range(self.n_from, self.n_to):
            plain = pd.read_csv(os.path.join(os.getcwd(), 'results', 'networks', 'conf_part'+str(j)+'_blockall_'+alg_sel+'_'+ct_sel+'_none_gene_list.csv'), header=0)
            for bl in list(self.conf_partitions[ct_sel][conf_sel].keys()):
                block_network = os.path.join(os.getcwd(), 'results', 'networks', 'rnd_part'+str(j)+'_block'+bl+'_'+alg_sel+'_'+ct_sel+'_'+conf_sel+'_gene_list.csv')
                inferred_networks = [plain, pd.read_csv(block_network, header=0)]
                for k in self.rep_k:
                    ji, state, s_int, s_un = mean_jaccard_index_at_k(k, inferred_networks, alg_sel)
                    rnd_results[j][bl].append(ji)
                JI_rnd = pd.DataFrame({'k': self.rep_k, 'mean JI': rnd_results[j][bl]}, 'partID':j, 'confounder':conf_sel,
                                        'method':alg_sel, 'cohort':ct_sel)
                networks_blocks_rnd[bl] = pd.concat(networks_blocks_rnd[bl], JI_rnd)

        for bl in block_ids:
            self.visualize_JIs(networks_blocks_conf[bl], networks_blocks_rnd[bl], bl)
            sign_conf_rnd = 0
            for k in self.rep_k:
                rnd = networks_blocks_rnd[bl][networks_blocks_rnd[bl]['k'] == k]
                conf = networks_blocks_conf[bl][networks_blocks_conf[bl]['k'] == k]
                wilcox_ = wilcoxon(conf['mean JI'], rnd['mean JI'], alternative='less', correction=True)
                if wilcox_.pvalue < 0.05:
                    sign_conf_rnd + = 1
            self._log(['One-sided wilcoxon test on ' + str(bl) + ' network and corresponding random network for ' +
                        str(sign_conf_rnd/len(self.rep_k)) + ' of all tested k.'])

        for bl0, bl1 in itt.product(block_ids, 2):
            net0 = networks_blocks_conf[bl0]
            net1 = networks_blocks_conf[bl1]
            sign_conf_conf = 0
            for k in self.rep_k:
                block0_conf = net0[net0['k'] == k]
                block1_conf = net1[net1['k'] == k]
                try:
                    wilcox = wilcoxon(block0_conf['mean JI'],block1_conf['mean JI'], alternative='less', correction=True)
                    if wilcox.pvalue < 0.05:
                        sign_conf_conf + = 1
                except:
                    assert all(x == 0.0 for x in (conf_0['mean JI']-conf_1['mean JI']))
                    continue
            
            self.wilcox_fractions = pd.concat([self.wilcox_fractions, pd.DataFrame({'cohort':ct_sel, 'confounder/variable':str(conf_sel), 
                            'method':str(alg_sel), 'bl0':str(bl0), 'bl1':str(bl1), 'frac':sign_conf_conf/len(self.rep_k)})])

    def visualize_mwus(ct_sel):
        plt.tight_layout()
        sns.set(font_scale=1.2)
        sns.set_style("darkgrid")
        palette={'age':"#2a60e1", 'ethnicity':"#e8e049", 'sex':"#69d447", 'WGCNA':'#15a5d9', 'CEMiTool':'#f9a62d', 'ARACNe-AP':'#d54f13'}
        frms = set(self.mwu_fractions['cohort'].values))
        h, axs_z = plt.subplots(frms,1, figsize=(14,5*frms), gridspec_kw=dict(width_ratios=[5], height_ratios=np.repeat(5, frms)))
        h.tight_layout(pad=5.0)

        for cohort in set(self.mwu_fractions['cohort'].values):
            cur = self.mwu_fractions[self.mwu_fractions['cohort'] == cohort]
            cur_var = pd.DataFrame({'method':cur['confounder/variable'], 'confounder/variable':cur['method'], 'frac':cur['frac']})
            cur = pd.concat([cur, cur_var])

            g =sns.swarmplot(ax=axs_z[i], x="method", y="frac",hue="confounder/variable",data=_cur, legend=True, s=15, palette=palette)
            g.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), ncol=1)
            ax = g.axes
            ax.yaxis.grid(True)
            ax.xaxis.grid(True)
            ax.axvline(2.5, ls='--')

            ax.set_title(cohort, loc='left', fontdict={'fontsize': 28,'verticalalignment': 'baseline','horizontalalignment': 'left'})
            g.set(xlabel='', ylabel = "fraction of k with \n significant p-value",title=cohort)
            ax.yaxis.label.set_size(18)
            g.tick_params(labelsize=16)
            plt.yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])

            i += 1

        fig = g.get_figure()
        h.savefig(os.path.join(cwd, 'mwu_comparison', str(cohort)+'_mwu_comparison.pdf'))

    def visualize_JIs(block_conf, block_rnd, bl):
        plt.figure()
        sns.set(font_scale=1.2)
        sns.set_style('whitegrid')
        block_conf['partition type'] = str(bl)
        block_rnd['partition type'] = 'random with size\nof ' + str(bl) + ' block'
        JI = pd.concat([block_conf, block_rnd])
        g = sns.lineplot(JI, x='k', y='mean JI', hue='partition type', errorbar='sd')
        g.set_title(alg_sel + ' ' + ct_sel + ' ' + conf_sel)
        plt.legend(bbox_to_anchor=(1.02, 0.9), loc='upper left', title='block compared with $G^{all}$', borderaxespad=0)
        #g.figure.savefig(os.path.join(os.getcwd(), 'res_plots', alg_sel+'_'+conf_sel+'_'+ct_sel+'_plain_compar.pdf'))

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
            path = os.path.join(cwd, 'results', 'networks', f'{mode}_part{part_nb}_block{block_id}_{alg_sel}_{ct_sel}_{conf_sel}_gene_list.csv')
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
                p, table = preprocessing._var_conf_chi(self.pheno_datasets[ct_sel], conf_sel, var, self.conf_dict)
                sign = ('**' if p < 0.01 else '*') if p < 0.05 else 'ns'
                self._log(['chi^2 test ' + str(conf_sel) +' and ' + var + ' significant: ' + str(p) + f'{sign}'])

    def _log(self, messages, mode='a'):
        if self.rank == 0:
            with open(self.logfile, mode) as f:
                for msg in messages:
                    f.write(msg)