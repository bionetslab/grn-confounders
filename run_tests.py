# -*- coding: utf-8 -*-

import test_suite
from test_suite import TestRunner
import argparse
import os
import itertools as itt
import traceback
from test_suite import Selectors
#from mpi4py import MPI

#comm = MPI.COMM_WORLD
#rank= comm.Get_rank()
#size = comm.Get_size()

def get_parser():
    """Return parser for command line argument processing."""
    parser = argparse.ArgumentParser('Assessing effect of sex, age, and ethnicity confounders on GRN and co-expression network inference.')
    parser.add_argument('-ct', required=True, nargs='+', choices=[str(sel) for sel in list(Selectors.CancerTypeSelector)])
    parser.add_argument('-conf', required=True, nargs='+', choices=[str(sel) for sel in list(Selectors.ConfounderSelector)])
    parser.add_argument('-alg', required=True, nargs='+', choices=[str(sel) for sel in list(Selectors.AlgorithmSelector)])
    parser.add_argument('-n', required=True, type=int)
    parser.add_argument('-m', required=True, type=int)
    parser.add_argument('-k', required=True, type=int)
    parser.add_argument('-mode', required=False, type=str)
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
    if args.mode == 'parallel':
        n_from = (int(args.n/size))*rank
        n_to = n_from + int(args.n/size)
        if rank == size-1:
            n_to = args.n
        m_from = (int(args.m/size))*rank
        m_to = m_from + int(args.m/size)
        if rank == size-1:
            m_to = args.m

        if rank == 0 and verbose:
            print('loading data ...')
        #test_runner = TestRunner(args.n, args.m, args.k)
        test_runner = TestRunner(n_from, n_to, m_from, m_to, args.k)
        if rank == 0 and verbose:
            print('running the tests ...')
        test_runner.run_on_cancer_types_confounders(args.ct, args.conf, args.alg, True)
    else:
        if verbose:
            print('loading data ...')
        test_runner = TestRunner(0, args.n, 0, args.m, args.k)
        if verbose:
            print('running the tests ...')
        test_runner.run_on_cancer_types_confounders(args.ct, args.conf, args.alg, True)
    return

if __name__ == '__main__':
    args = get_parser().parse_args()
    run_tests(args)
    
