# -*- coding: utf-8 -*-
import test_suite
from test_suite.TestRunner import TestRunner
import argparse
import os
import itertools as itt
import traceback
from test_suite import Selectors
import yaml

def get_parser():
    """Return parser for command line argument processing."""
    parser = argparse.ArgumentParser('Assessing effect of sex, age, and ethnicity confounders on GRN and co-expression network inference.')
    data_parser = parser.add_subparsers(dest='input', required=True, help='Use command line (cmd) or config files (config) to specify input parameters.')

    config = data_parser.add_parser('config', help='Use config file data.yml, fields.yml, params.yml to specify input parameters.')

    custom = data_parser.add_parser('cmd', help='Use command line to specify input parameters.') #'Use own data. Place data in data directory and specify filenames in -ged and -pt. Order of files in -ged must match order of files in -pt, i.e. the ct, ged, and pt at index 0 correspond to each other.')
    custom.add_argument('-tcga', action='store_true', help='run tests additionally on combined datasets.')
    
    custom.add_argument('-ct', required=True, nargs='+', help='Define cohort names for the filenames in ged and pt.')
    custom.add_argument('-ged', required=False, nargs='*', help='Specify filename of gene expression data.')
    custom.add_argument('-pt', required=False, nargs='*', help='Specify filename of pheno type data.')
    custom.add_argument('-sep', required=False, nargs='*', help='Separator in pt and ged files. Default is \',\'.')

    custom.add_argument('-conf', required=False, nargs='*')
    custom.add_argument('-block_types', required=False, nargs='*', choices=[str(sel) for sel in list(Selectors.BlockType)])
    custom.add_argument('-roles', required=False, nargs='*')

    custom.add_argument('-alg', required=True, nargs='+', choices=[str(sel) for sel in list(Selectors.AlgorithmSelector)])
    custom.add_argument('-N_from', required=False, type=int, nargs='?', const=0, default=0)
    custom.add_argument('-M_from', required=False, type=int, nargs='?', const=0, default=0)
    custom.add_argument('-N_to', required=False, type=int, nargs='?', const=20, default=20)
    custom.add_argument('-M_to', required=False, type=int, nargs='?', const=20, default=20)
    custom.add_argument('-k', required=False, type=int, nargs='?', const=5000, default=5000)
    custom.add_argument('-combine', action='store_true', help='run tests additionally on combined datasets.')
    custom.add_argument('-par', action='store_true', help='If set, run tests in parallel, else, run tests sequentially.')
    custom.add_argument('-g_all', action='store_true', help='If set, run method on entire dataset and infer g_all.')

    return parser

def dump_config(args):
    """instantiates testRunner object and passes the user arguments.

    Parameters
    ----------
    args: argparse.Namespace
        Namespace object populated with user command line arguments.
    """
    assert len(args.block_types) == len(args.conf) == len(args.roles) if args.conf is not None else True
    assert len(args.conf) > 0 or args.g_all, 'If -g_all flag is not set, specify at least one confounder.'
    assert len(args.ct) == len(args.ged) == len(args.pt) == len(args.sep) if not args.tcga else True

    data = dict()
    for ct in args.ct:
        if args.tcga:
            assert ct in [str(el) for el in list(Selectors.TCGACancerTypeSelector)]
            data[ct] = {'tcga': True, 'ged': None, 'pt': None, 'sep': None}
        else:
            g, p, s = args.ged, args.pt, args.sep
            data[ct] = {'tcga': False, 'ged': g, 'pt': p, 'sep': s}

    params = {'algorithms': args.alg, 'N_from': args.N_from, 'N_to': args.N_to, 'M_from': args.M_from, 'M_to': args.M_to,
                'k': args.k, 'combine': args.combine, 'par': args.par, 'g_all': args.g_all}

    fields = dict()
    for (conf_sel, conf_role, conf_type) in zip(args.conf, args.roles, args.block_types):
        fields[conf_sel] = {'role': conf_role, 'type': conf_type}

    cwd = os.getcwd()
    with open(os.path.join(cwd, 'config', 'data_cmd.yml'), 'w') as f:
        data = yaml.dump(data, f)
    with open(os.path.join(cwd, 'config', 'fields_cmd.yml'), 'w') as f:
        data = yaml.dump(fields, f)
    with open(os.path.join(cwd, 'config', 'params_cmd.yml'), 'w') as f:
        data = yaml.dump(params, f)

def parse_config(data_p, fields_p, params_p):
    # read parameters from .yml file
    cwd = os.getcwd()
    config_path = os.path.join(cwd, 'config')

    data = dict()
    with open(os.path.join(config_path, data_p)) as f:
        _data = yaml.load_all(f, Loader=yaml.FullLoader)
        for doc in _data:
            data = doc
            break
        for key in data.keys():
            if data[key]['tcga']:
                data[key] = {'ged': f'TCGA-{key}.htseq_fpkm.tsv', 'pt': f'TCGA-{key}.GDC_phenotype.tsv', 'sep': ','}
            else:
                assert data[key]['ged'] is not None and data[key]['ged'] is not None and data[key]['sep'], 'Specify ged (gene expression data) file name, pt (pheno type) file name, and sep (separators) for each cohort if tcga option is set to False.'

    params = {'algorithm':'', 'N_from': '', 'N_to': '', 'M_from': '', 'M_to': '', 'k': '', 'combine': '', 'par': '', 'g_all': ''}
    with open(os.path.join(config_path, params_p)) as f:
        _params = yaml.load_all(f, Loader=yaml.FullLoader)
        for doc in _params:
            assert len(doc) == len(params.keys())
            params = doc
            break

    fields = dict()
    with open(os.path.join(config_path, fields_p)) as f:
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
    
    return data, fields, params

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
                            g_all=params['g_all'], rank=rank)
    test_runner.run_on_cancer_types_confounders()

if __name__ == '__main__':
    args = get_parser().parse_args()
    if args.input == 'cmd':
        dump_config(args)
        data_p, fields_p, params_p = 'data_cmd.yml', 'fields_cmd.yml', 'params_cmd.yml'
    elif args.input == 'config':
        data_p, fields_p, params_p = 'data.yml', 'fields.yml', 'params.yml'
    data, fields, params = parse_config(data_p, fields_p, params_p)
    run_tests(data, fields, params)