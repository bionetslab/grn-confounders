import test_suite
from test_suite import TestRunner
import argparse
import os
import itertools as itt
import traceback
#from test_suite import Selectors

def run_tests(verbose=True):
    #try:
    if verbose:
        print('loading data ...')
    test_runner = TestRunner(1)
    #print(test_runner.conf_partitions[Selectors.CancerTypeSelector.BLCA].keys())
    if verbose:
        print('running the tests ...')
    test_runner.run_on_all_cancer_types_confounders_partitions(verbose=True)
    if verbose:
        print('saving the results')
    #test_runner.save_results()
    return 0
    #except Exception:
        #traceback.print_exc()
        #return 1

if __name__ == '__main__':
    cwd = os.path.join(os.path.dirname(__file__))
    os.chdir(cwd)
    print(os.getcwd())
    run_tests()
    
