# -*- coding: utf-8 -*-
from confinspect import TestRunner
import InputHandler
import os
import argparse

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

if __name__ == '__main__':
    #InputHandler.setup_directories()
    #args = get_parser().parse_args()
    #if args.input == 'cmd':
        #data_p, fields_p, params_p = InputHandler.dump_config(args)
    #elif args.input == 'config':
    data_p, fields_p, params_p = 'data.yml', 'fields.yml', 'params.yml'
    data, fields, params = InputHandler.parse_config(data_p, fields_p, params_p)
    run_tests(data, fields, params)
