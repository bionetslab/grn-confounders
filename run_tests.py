# -*- coding: utf-8 -*-
import test_suite
from test_suite import TestRunner
import argparse
import os
import itertools as itt
import traceback
from test_suite import Selectors

def get_parser():
    """Return parser for command line argument processing."""
    parser = argparse.ArgumentParser('Assessing effect of sex, age, and ethnicity confounders on GRN and co-expression network inference.')
    data_parser = parser.add_subparsers(dest='data', required=True, help='Use TCGA data or custom data.')

    parser.add_argument('-conf', required=True, nargs='+')
    parser.add_argument('-block_types', required=True, nargs='+', choices=[str(sel) for sel in list(Selectors.BlockType)])
    parser.add_argument('-alg', required=True, nargs='+', choices=[str(sel) for sel in list(Selectors.AlgorithmSelector)])
    parser.add_argument('-k', required=True, type=int)
    parser.add_argument('-combine', action='store_true', help='run tests additionally on combined datasets.')
    parser.add_argument('par', action='store_true', help='If set, run tests in parallel, else, run tests sequentially.')
    parser.add_argument('g_all', action='store_true', help='If set, run method on entire dataset and infer g_all.')
    parser.add_argument('-sep', required=False, const='\t', type=str, help='Separator in pt and ged files. Default is \'\\t\'.')
    parser.add_argument('-chi', nargs='+')

    tcga = data_parser.add_parser('tcga', help='Use TCGA data. Specify study abbreviation of TCGA studies as cancer types in -ct.')
    tcga.add_argument('-ct', required=True, nargs='+')

    custom = data_parser.add_parser('custom', help='Use own data. Place data in data directory and specify filenames in -ged and -pt.\
         Order of files in -ged must match order of files in -pt, i.e. the ct, ged, and pt at index 0 correspond to each other.')
    custom.add_argument('-ct', required=True, nargs='+', help='Define cohort names for the filenames in ged and pt.')
    custom.add_argument('-ged', required=True, nargs='+', help='Specify filename of gene expression data.')
    custom.add_argument('-pt', required=True, nargs='+', help='Specify filename of pheno type data.')

    return parser

def run_tests(args):
    """instantiates testRunner object and passes the user arguments.

    Parameters
    ----------
    args: argparse.Namespace
        Namespace object populated with user command line arguments.
    """
    assert args.N_from < args.N_to, 'N_from must be smaller than N_to'
    assert args.M_from < args.M_to, 'M_from must be smaller than M_to'
    assert len(args.block_types) == len(args.conf)
    assert len(args.ct) == len(args.ged) == len(args.pt)
    assert all([el in args.conf for el in args.chi])

    # create data_dict with paths to gene expression data and phenotype data
    data_dict = {c:{} for c in args.ct}
    if args.data == 'tcga':
        sep = '\t'
        for i in range(len(args.ct)):
            data_path[args.ct[i]] = {'ged': 'TCGA-{ct[i]}.htseq_fpkm.csv', 'pt': 'TCGA-{ct[i]}.GDC_phenotype.csv'}
    elif args.data == 'custom':
        sep = args.sep
        for i in range(len(args.ct)):
            data_dict[args.ct[i]] = {'ged': args.ged[i], 'pt': args.pt[i]}

    # create conf_dict with confounders and corresponding block types (needed for generation of confounder-induced partition)
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
    test_runner = TestRunner(data_dict, args.conf, args.alg, n_from, n_to, m_from, m_to, args.k, combine=args.combine, sep=sep, 
                            g_all=args.g_all, chi_variables=args.chi, rank=rank)
    test_runner.run_on_cancer_types_confounders(args.combine, True)
    return

if __name__ == '__main__':
    args = get_parser().parse_args()
    run_tests(args)
    
