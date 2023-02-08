# -*- coding: utf-8 -*-
from confinspect import TestRunner
import InputHandler
import os
import argparse
import yaml

def get_parser():
    """Return parser for command line argument processing."""
    parser = argparse.ArgumentParser('Assessing effect of sex, age, and ethnicity confounders on GRN and co-expression network inference.')
    data_parser = parser.add_subparsers(dest='input', required=True, help='Use command line (cmd) or config files (config) to specify input parameters.')

    config = data_parser.add_parser('config', help='Use config file data.yml, fields.yml, params.yml to specify input parameters.')

    custom = data_parser.add_parser('cmd', help='Use command line to specify input parameters.')
    custom.add_argument('-tcga', action='store_true', help='run tests additionally on combined datasets.')

    custom.add_argument('-ct', required=True, nargs='+', help='Define cohort names for the filenames in ged and pt.')
    custom.add_argument('-ged', required=False, nargs='*', help='Specify filename of gene expression data.')
    custom.add_argument('-pt', required=False, nargs='*', help='Specify filename of pheno type data.')
    custom.add_argument('-sep', required=False, nargs='*', help='Separator in pt and ged files. Default is \',\'.')

    custom.add_argument('-conf', required=False, nargs='*', const=[], default=[])
    custom.add_argument('-block_types', required=False, nargs='*', choices=[str(sel) for sel in list(Selectors.BlockType)], const=[], default=[])
    custom.add_argument('-roles', required=False, nargs='*',  choices=[str(sel) for sel in list(Selectors.Role)], const=[], default=[])

    custom.add_argument('-alg', required=True, nargs='+', choices=[str(sel) for sel in list(Selectors.AlgorithmSelector)])
    custom.add_argument('-N_from', required=False, type=int, nargs='?', const=0, default=0)
    custom.add_argument('-M_from', required=False, type=int, nargs='?', const=0, default=0)
    custom.add_argument('-N_to', required=False, type=int, nargs='?', const=100, default=100)
    custom.add_argument('-M_to', required=False, type=int, nargs='?', const=10, default=10)
    custom.add_argument('-k_max', required=False, type=int, nargs='?', const=5000, default=5000)
    custom.add_argument('-tissue_type_field', required=False, type=int, nargs='?', const=None, default=None)
    custom.add_argument('-tissue_type', required=False, type=int, nargs='?', const=None, default=None)
    custom.add_argument('-combine', action='store_true', help='run tests additionally on combined datasets.')
    custom.add_argument('-par', action='store_true', help='If set, run tests in parallel, else, run tests sequentially.')
    custom.add_argument('-g_all', action='store_true', help='If set, run method on entire dataset and infer g_all.')
    custom.add_argument('-save_networks', action='store_true', help='If set, save generated networks. Caution: memory intensive!')
    custom.add_argument('-log', required=False, type=str, const='log.txt', default='log.txt')
    return parser

def run_tests(data, fields, params):
    """Instantiates TestRunner object and passes the user arguments.
    Parameters
    ----------
    data : dict
        Dictionary populated with user input about data.
    fields : dict
        Dictionary populated with user input about confounder and variable fields.
    params : dict
        Dictionary populated with user input about parameters.
    """
    InputHandler.verify_input(data, params, fields)
    # prepare parallel vs. sequential
    if params['par']:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        rank= comm.Get_rank()
        size = comm.Get_size()
    else:
        rank = 0
        size = 1

    # set starting and ending indices of workers
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
    test_runner = TestRunner.TestRunner(data, fields, params, rank=rank)
    test_runner.induce_partitions()
    test_runner.run_all()

def dump_config(args):
    """Dump cmd line input arguments to config files for documentation and reusability purposes. Reuse code for reading config file
    afterwards.
    Parameters
    ----------
    args: argparse.Namespace
        Namespace object populated with user command line arguments.
    """
    assert len(args.block_types) == len(args.conf) == len(args.roles) if args.conf is not None else True
    assert len(args.ct) == len(args.ged) == len(args.pt) == len(args.sep) if not args.tcga else True
    
    cwd = os.getcwd()
    data = dict()
    for ct in args.ct:
        if args.tcga:
            data[ct] = {'tcga': True, 'ged': None, 'pt': None, 'sep': None, 'tissue_type_field': None, 'tissue_type': None}
        else:
            g, p, s, tf, tt = args.ged, args.pt, args.sep, args.tissue_type_field, args.tissue_type
            data[ct] = {'tcga': False, 'ged': g, 'pt': p, 'sep': s, 'tissue_type_field': tf, 'tissue_type': tt}

    params = {'algorithms': args.alg, 'N_from': args.N_from, 'N_to': args.N_to, 'M_from': args.M_from, 'M_to': args.M_to,
                'k_max': args.k_max, 'combine': args.combine, 'par': args.par, 'g_all': args.g_all, 'save_networks': args.save_networks, 'logfile': args.log}

    fields = dict()
    for (conf_sel, conf_role, conf_type) in zip(args.conf, args.roles, args.block_types):
        fields[conf_sel] = {'role': conf_role, 'type': conf_type}

    if not os.path.exists(config_path):
        os.mkdir(config_path)
    with open(os.path.join(cwd, 'config', 'data_cmd.yml'), 'w') as f:
        data = yaml.dump(data, f)
    with open(os.path.join(cwd, 'config', 'fields_cmd.yml'), 'w') as f:
        data = yaml.dump(fields, f)
    with open(os.path.join(cwd, 'config', 'params_cmd.yml'), 'w') as f:
        data = yaml.dump(params, f)

    return 'data_cmd.yml', 'fields_cmd.yml', 'params_cmd.yml'

def parse_config(data_p, fields_p, params_p):
    cwd = os.getcwd()
    config_path = os.path.join(cwd, 'config')
    assert os.path.exists(config_path) and os.path.exists(os.path.join(config_path, data_p)) and os.path.exists(os.path.join(config_path, fields_p)) and os.path.exists(os.path.join(config_path, params_p)), 'Put data.yml, fields.yml, and params.yml in config/'
    
    data = dict()
    with open(os.path.join(config_path, data_p)) as f:
        _data = yaml.load_all(f, Loader=yaml.FullLoader)
        for doc in _data:
            data = doc
            break
    with open(os.path.join(config_path, params_p)) as f:
        _params = yaml.load_all(f, Loader=yaml.FullLoader)
        for doc in _params:
            params = doc
            break
    fields = dict()
    with open(os.path.join(config_path, fields_p)) as f:
        _fields = yaml.load_all(f, Loader=yaml.FullLoader)
        for doc in _fields:
            fields = doc
            break
    return data, fields, params

if __name__ == '__main__':
    #InputHandler.setup_directories()
    #args = get_parser().parse_args()
    #if args.input == 'cmd':
        #data_p, fields_p, params_p = dump_config(args)
    #elif args.input == 'config':
    data_p, fields_p, params_p = 'data.yml', 'fields.yml', 'params.yml'
    data, fields, params = parse_config(data_p, fields_p, params_p)
    run_tests(data, fields, params)
