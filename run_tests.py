# -*- coding: utf-8 -*-

import test_suite
from test_suite import TestRunner
import argparse
import os
import itertools as itt
import traceback
import test_suite

def get_parser():
    """TODO"""
    parser = argparse.ArgumentParser('Assessing effect of sex, age, and ethnicity confounders on GRN and co-expression network inference.')
    parser.add_argument(dest=alg_list, required=True, nargs='+', choices=list(Selectors.AlgorithmSelectors))
    parser.add_argument(dest=conf_list, required=True, nargs='+', choices=list(Selectors.ConfounderSelectors))
    parser.add_argument(dest=alg_list, required=True, nargs='+', choices=list(Selectors.AlgorithmSelectors))
    parser.add_argument(dest=n, required=True)
    parser.add_argument(dest=k, required=True)
    return parser

def run_tests(n, verbose=True):
    """runs the tests.

    Parameters
    ----------
    cancer_types : list
        List of strings that specify the cohorts that should be investigated.

    confounders : list
        List of strings that specify the confounders that should be investigated.

    verbose : bool
        Print progress to stdout.
    """
    args = get_parser().parse_args()

    if verbose:
        print('loading data ...')
    test_runner = TestRunner(args.n, args.k)
    if verbose:
        print('running the tests ...')
    test_runner.run_on_cancer_types_confounders(args.alg_list, args.conf_list, args.alg_list, True)
    return 0

if __name__ == '__main__':
    cwd = os.path.join(os.path.dirname(__file__))
    os.chdir(cwd) # on HPC, this line needs to be commented, as the batch_script cd-s into the right directory

    run_tests(3)
    
