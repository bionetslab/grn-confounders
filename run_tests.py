# -*- coding: utf-8 -*-
from confinspect import TestRunner
import InputHandler
import os
import yaml

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
    InputHandler.setup_directories()
    data_p, fields_p, params_p = 'data.yml', 'fields.yml', 'params.yml'
    data, fields, params = parse_config(data_p, fields_p, params_p)
    print(params)
    print(fields)
    print(data)
    run_tests(data, fields, params)
