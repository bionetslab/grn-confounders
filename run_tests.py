# -*- coding: utf-8 -*-

import test_suite
from test_suite import TestRunner
import argparse
import os
import itertools as itt
import traceback
from test_suite import Selectors
"""
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank= comm.Get_rank()
size = comm.Get_size()
"""
def get_parser():
    """Return parser for command line argument processing.
    E.g.: python run_tests.py -ct BLCA -conf age -alg ARACNE -k 500  par -N_from 0 -N_to 1000 -M_from 0 -M_to 100 """
    parser = argparse.ArgumentParser('Assessing effect of sex, age, and ethnicity confounders on GRN and co-expression network inference.')
    mode_parser = parser.add_subparsers(dest='mode', required=True, help='run test sequentially or in parallel.')

    parser.add_argument('-ct', required=True, nargs='+', choices=[str(sel) for sel in list(Selectors.CancerTypeSelector)])
    parser.add_argument('-conf', required=True, nargs='+', choices=[str(sel) for sel in list(Selectors.ConfounderSelector)])
    parser.add_argument('-alg', required=True, nargs='+', choices=[str(sel) for sel in list(Selectors.AlgorithmSelector)])
    parser.add_argument('-k', required=True, type=int)
    
    seq = mode_parser.add_parser('seq', help='run tests sequentially.')
    seq.add_argument('-n', required=True, type=int)
    seq.add_argument('-m', required=True, type=int)

    par = mode_parser.add_parser('par', help='run tests in parallel.')
    par.add_argument('-N_from', required=True, type=int) # 0 to run full batch of partitions in parallel
    par.add_argument('-N_to', required=True, type=int)
    par.add_argument('-M_from', required=True, type=int) # 0 to run full batch of partitions in parallel
    par.add_argument('-M_to', required=True, type=int)

    return parser

def run_tests(args, verbose=True):
    """instantiates testRunner object and passes the user arguments.

    Parameters
    ----------
    args: argparse.Namespace
        Namespace object populated with user command line arguments.

    verbose : bool
        Print progress to stdout.
    """
    if args.mode == 'par':
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

        if rank == 0 and verbose:
            print('loading data ...')
        test_runner = TestRunner(n_from, n_to, m_from, m_to, args.k, rank)
        if rank == 0 and verbose:
            print('running the tests ...')
        test_runner.run_on_cancer_types_confounders(args.ct, args.conf, args.alg, True)
    else:
        if verbose:
            print('loading data ...')
        test_runner = TestRunner(0, args.n, 0, args.m, args.k, 0)
        if verbose:
            print('running the tests ...')
        test_runner.run_on_cancer_types_confounders(args.ct, args.conf, args.alg, True)
    return

if __name__ == '__main__':
    args = get_parser().parse_args()
    run_tests(args)
    
