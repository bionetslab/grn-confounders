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
    parser.add_argument('-ct', required=True, nargs='+', choices=[str(sel) for sel in list(Selectors.CancerTypeSelector)])
    parser.add_argument('-conf', required=True, nargs='+', choices=[str(sel) for sel in list(Selectors.ConfounderSelector)])
    parser.add_argument('-alg', required=True, nargs='+', choices=[str(sel) for sel in list(Selectors.AlgorithmSelector)])
    parser.add_argument('-n', required=True, type=int)
    parser.add_argument('-k', required=True, type=int)
    return parser

def run_tests(args, verbose=True):
    """runs the tests.

    Parameters
    ----------
    args: Namespace
        Namespace object populated with user command line arguments.

    verbose : bool
        Print progress to stdout.
    """

    if verbose:
        print('loading data ...')
    test_runner = TestRunner(args.n, args.k)
    if verbose:
        print('running the tests ...')
    test_runner.run_on_cancer_types_confounders(args.ct, args.conf, args.alg, True)
    return 0

if __name__ == '__main__':
    args = get_parser().parse_args()
    cwd = os.getcwd()
    #os.chdir(cwd) # on HPC, this line needs to be commented, as the batch_script cd-s into the right directory
    run_tests(args)
    
