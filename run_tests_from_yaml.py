# -*- coding: utf-8 -*-
import test_suite
from test_suite.TestRunner import TestRunner
import argparse
import os
import itertools as itt
import traceback
from test_suite import Selectors
import yaml
import pprint

def run_tests(data, fields, params):
    """instantiates testRunner object and passes the user arguments.

    Parameters
    ----------
    data : dict
        Dictionary populated with user input about data.
    fields : dict
        Dictionary populated with user input about confounder and variable fields.
    params : dict
        Dictionary populated with user input about parameters.
    """
    assert params['N_from'] <= params['N_to'], 'N_from must be smaller than or equal to N_to'
    assert params['M_from'] <= params['M_to'], 'M_from must be smaller than or equal to M_to'
    assert params['g_all'] or len(list(data.keys())) > 0, 'If no confounders specified, -g_all flag must be set to infer network from entire data.'
    chi_variables = [key for key in fields.keys() if fields[key]['role'] == 'variable']

    # prepare parallel vs. sequential
    if params['par']:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        rank= comm.Get_rank()
        size = comm.Get_size()
    else:
        rank = 0
        size = 1

    # set starting and ending points of workers
    batch_size_n = int((params['N_to'] - params['N_from'])/size)
    n_from = params['N_from'] + (batch_size_n)*rank 
    n_to = n_from + batch_size_n
    if rank == size-1:
        n_to = params['N_to']
    batch_size_m = int((params['M_to'] - params['M_from'])/size)
    m_from = params['M_from'] + batch_size_m*rank
    m_to = m_from + batch_size_m
    if rank == size-1:
        m_to = params['M_to']

    # start tests
    test_runner = TestRunner(data, fields, params['algorithms'], n_from, n_to, m_from, m_to, params['k'], combine=params['combine'],
                            g_all=params['g_all'], chi_variables=chi_variables, rank=rank)
    test_runner.run_on_cancer_types_confounders()
    return

if __name__ == '__main__':
    # read parameters from .yml file
    cwd = os.getcwd()
    config_path = os.path.join(cwd, 'config')

    data = dict()
    with open(os.path.join(config_path, 'data.yml')) as f:
        _data = yaml.load_all(f, Loader=yaml.FullLoader)
        for doc in _data:
            data = doc
            break
        for key in data.keys():
            if data[key]['tcga']:
                data[key] = {'ged': f'TCGA-{key}.htseq_fpkm.tsv', 'pt': f'TCGA-{key}.GDC_phenotype.tsv', 'sep': ','}
            else:
                assert data[key]['ged'] is not None, 'Specify ged (gene expression data) file name and pt (pheno type) file name if tcga option is set to False.'

    params = {'algorithm':'', 'N_from': '', 'N_to': '', 'M_from': '', 'M_to': '', 'k': '', 'combine': '', 'par': '', 'g_all': ''}
    with open(os.path.join(config_path, 'params.yml')) as f:
        _params = yaml.load_all(f, Loader=yaml.FullLoader)
        for doc in _params:
            assert len(doc) == len(params.keys())
            params = doc
            break

    fields = dict()
    with open(os.path.join(config_path, 'fields.yml')) as f:
        _fields = yaml.load_all(f, Loader=yaml.FullLoader)
        for doc in _fields:
            fields = doc
            break
        for key in fields.keys():
            if key is None:
                assert params['g_all'] == True, 'If no fields are specified, g_all must be set to True. Otherwise, nothing to do here.'

        # convert types of confounders and variables into Selector identifier
        for key in fields.keys():
            fields[key]['type'] = Selectors.BlockType(fields[key]['type'])
        # if user wishes to compare with G_all, add 'NONE' to confounders
        if params['g_all']:
            fields[Selectors.ConfounderSelector.NONE] = {'role': 'confounder', 'type': Selectors.BlockType.ALL}

    run_tests(data, fields, params)
    
