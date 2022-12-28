# -*- coding: utf-8 -*-
import test_suite
from test_suite.TestRunner import TestRunner
import argparse
import os
import itertools as itt
import traceback
from test_suite import Selectors

def get_parser():
    """Return parser for command line argument processing."""
    parser = argparse.ArgumentParser('Assessing effect of sex, age, and ethnicity confounders on GRN and co-expression network inference.')
    data_parser = parser.add_subparsers(dest='data', required=True, help='Use TCGA data or custom data.')

    tcga = data_parser.add_parser('tcga', help='Use TCGA data. Specify study abbreviation of TCGA studies as cancer types in -ct.')
    tcga.add_argument('-ct', required=True, nargs='+', choices=[str(sel) for sel in list(Selectors.TCGACancerTypeSelector)])

    tcga.add_argument('-N_from', required=False, type=int, nargs='?', const=0, default=0)
    tcga.add_argument('-M_from', required=False, type=int, nargs='?', const=0, default=0)
    tcga.add_argument('-N_to', required=False, type=int, nargs='?', const=20, default=20)
    tcga.add_argument('-M_to', required=False, type=int, nargs='?', const=20, default=20)

    tcga.add_argument('-conf', required=False, nargs='*')
    tcga.add_argument('-block_types', required=True, nargs='+', choices=[str(sel) for sel in list(Selectors.BlockType)])
    tcga.add_argument('-alg', required=True, nargs='+', choices=[str(sel) for sel in list(Selectors.AlgorithmSelector)])
    tcga.add_argument('-k', required=False, type=int, nargs='?', const=5000, default=5000)
    tcga.add_argument('-combine', action='store_true', help='run tests additionally on combined datasets.')
    tcga.add_argument('-par', action='store_true', help='If set, run tests in parallel, else, run tests sequentially.')
    tcga.add_argument('-g_all', action='store_true', help='If set, run method on entire dataset and infer g_all.')
    tcga.add_argument('-sep', required=False, nargs='?', const='\t', default='\t', type=str, help='Separator in pt and ged files. Default is \'\\t\'.')
    tcga.add_argument('-chi', nargs='+')

    custom = data_parser.add_parser('custom', help='Use own data. Place data in data directory and specify filenames in -ged and -pt.\
         Order of files in -ged must match order of files in -pt, i.e. the ct, ged, and pt at index 0 correspond to each other.')
    custom.add_argument('-ct', required=True, nargs='+', help='Define cohort names for the filenames in ged and pt.')
    custom.add_argument('-ged', required=True, nargs='+', help='Specify filename of gene expression data.')
    custom.add_argument('-pt', required=True, nargs='+', help='Specify filename of pheno type data.')

    custom.add_argument('-N_from', required=False, type=int, nargs='?', const=0, default=0)
    custom.add_argument('-M_from', required=False, type=int, nargs='?', const=0, default=0)
    custom.add_argument('-N_to', required=False, type=int, nargs='?', const=20, default=20)
    custom.add_argument('-M_to', required=False, type=int, nargs='?', const=20, default=20)

    custom.add_argument('-conf', required=False, nargs='*')
    custom.add_argument('-block_types', required=False, nargs='+', choices=[str(sel) for sel in list(Selectors.BlockType)])
    custom.add_argument('-alg', required=True, nargs='+', choices=[str(sel) for sel in list(Selectors.AlgorithmSelector)])
    custom.add_argument('-k', required=False, type=int, nargs='?', const=5000, default=5000)
    custom.add_argument('-combine', action='store_true', help='run tests additionally on combined datasets.')
    custom.add_argument('-par', action='store_true', help='If set, run tests in parallel, else, run tests sequentially.')
    custom.add_argument('-g_all', action='store_true', help='If set, run method on entire dataset and infer g_all.')
    custom.add_argument('-sep', required=False, nargs='?', const='\t', default='\t', type=str, help='Separator in pt and ged files. Default is \'\\t\'.')
    custom.add_argument('-chi', nargs='+')


    return parser

def run_tests(args):
    """instantiates testRunner object and passes the user arguments.

    Parameters
    ----------
    args: argparse.Namespace
        Namespace object populated with user command line arguments.
    """
    assert args.N_from <= args.N_to, 'N_from must be smaller than or equal to N_to'
    assert args.M_from <= args.M_to, 'M_from must be smaller than or equal to M_to'
    assert len(args.block_types) == len(args.conf) if args.conf is not None else True
    assert all([el in args.conf for el in args.chi]) if args.chi is not None else True
    assert args.g_all or len(list(args.conf)) > 0, 'If no confounders specified, -g_all flag must be set to infer network from entire data.'
    
    # create data_dict with paths to gene expression data and phenotype data
    data_dict = {c:{} for c in args.ct}
    types = list(args.ct)
    if args.data == 'tcga':
        sep = '\t'
        for i in range(len(types)):
            data_dict[types[i]] = {'ged': f'TCGA-{types[i]}.htseq_fpkm.tsv', 'pt': f'TCGA-{types[i]}.GDC_phenotype.tsv'}
    elif args.data == 'custom':
        ged = list(args.ged)
        pt = list(args.pt)
        assert len(types) == len(ged) == len(pt)
        sep = args.sep
        for i in range(len(types)):
            data_dict[types[i]] = {'ged': ged[i], 'pt': pt[i]}

    # create conf_dict with confounders and corresponding block types (needed for generation of confounder-induced partition)
    args.conf = {} if args.conf is None else args.conf
    conf_dict = {c:'' for c in args.conf}
    for i in range(len(args.conf)):
        conf_dict[args.conf[i]] = Selectors.BlockType(args.block_types[i])
    # if user wishes to compare with G_all, add 'NONE' to confounders
    if args.g_all:
        conf_dict[Selectors.ConfounderSelector.NONE] = Selectors.BlockType.ALL

    # prepare parallel vs. sequential
    if args.par:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        rank= comm.Get_rank()
        size = comm.Get_size()
    else:
        rank = 0
        size = 1

    # set starting and ending points of workers
    batch_size_n = int((args.N_to - args.N_from)/size)
    n_from = args.N_from + (batch_size_n)*rank 
    n_to = n_from + batch_size_n
    if rank == size-1:
        n_to = args.N_to
    batch_size_m = int((args.M_to - args.M_from)/size)
    m_from = args.M_from + batch_size_m*rank
    m_to = m_from + batch_size_m
    if rank == size-1:
        m_to = args.M_to

    # start tests
    test_runner = TestRunner(data_dict, conf_dict, args.alg, n_from, n_to, m_from, m_to, args.k, combine=args.combine, sep=sep, 
                            g_all=args.g_all, chi_variables=args.chi, rank=rank)
    test_runner.run_on_cancer_types_confounders()
    return

if __name__ == '__main__':
    args = get_parser().parse_args()
    run_tests(args)
    
