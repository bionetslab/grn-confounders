import argparse
import yaml
import os
from confinspect import Selectors

def setup_directories():
    cwd = os.getcwd()
    dirs = [
    os.path.join(cwd, 'config'),
    os.path.join(cwd, 'data'),
    os.path.join(cwd, 'results'),
    os.path.join(cwd, 'results', 'JI'),
    os.path.join(cwd, 'results', 'networks'),
    os.path.join(cwd, 'data'),
    os.path.join(cwd, 'plots'),
    os.path.join(cwd, 'temp'),
    os.path.join(cwd, 'partitions')
    ]
    for dir_ in dirs:
        if not os.path.exists(dir_):
            os.mkdir(dir_)
    algrithms = os.path.join(cwd, 'algorithms')
    if not os.path.exists(algrithms):
        print('WARNING: directory with caller scripts for network inference algorithms not present. Download the directory from TODO if you want to \
            to use the predefined NetworkInferenceWrappers to run any of the algorithms ARACNe-AP, CEMiTool, GENIE3, GRNBoost2, or WGCNA. \
            Follow the installation instructions given in the README at TODO. For further instructions, see the README at TODO.')
        print('INFO: the code expects the user to implement a custom NetworkInferenceWrapper. Inherit from the abstract NetworkInferenceWrapper \
            of this package and hand the wrapper to the TestRunner object manually by calling add_custom_algorithm.')

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

def default_missing_data(data):
    """Check data dict. E.g. data_dict = {'HNSC': {'tcga': False, 'ged': 'ged_filename.csv', 'pt': 'pt_filename.csv', 'sep': ',', 'tissue_type_field': 'sample_type.sample', 'tissue_type': 'Primary Tumor'}}. Further documentation on the data_dict can be found at TODO.
    Parameters
    ----------
    data : dict
        Dictionary populated with user input about data.
    Return
    ----------
    data : dict
        Checked data dict.
    """
    for ct_dict in data.values():
        ct_dict['tcga'] = False if 'tcga' not in ct_dict.keys() else ct_dict['tcga']
        if 'sep' not in ct_dict.keys():
            print('Setting sep separator to \',\'...')
            ct_dict['sep'] = ','
        if 'tissue_type_field' not in ct_dict.keys():
            ct_dict['tissue_type_field'] = None
            ct_dict['tissue_type'] = None
        if 'tissue_type' not in ct_dict.keys():
            ct_dict['tissue_type_field'] = None
            ct_dict['tissue_type'] = None
    return data

def default_missing_params(params):
    """Check params dict. E.g. params_dict = {'algorithms':'', 'N_from': '', 'N_to': '', 'M_from': '', 'M_to': '', 'k_max': '', 'combine': '', 'par': '', 'g_all': '', 'save_networks': '', 'logfile': ''}. Further documentation on the params_dict can be found at TODO.
    Add default values N_from, M_from = 0, N_to = N_from + 100, M_to = M_from + 10, k_max = 5000, combine, par, g_all, save_networks = False if missing in params dict.
    Parameters
    ----------
    params : dict
        Dictionary populated with user input about parameters.
    Return
    ----------
    params : dict
        Checked params dict.
    """
    if 'algorithms' not in params.keys():
        print('No algorithm specified in params[\'algorithms\']. Please add custom algorithm after instantiating TestRunner by calling TODO before calling induce_partitions().')
        params['algorithms'] = []
    params['N_from'] = 0 if 'N_from' not in params.keys() else params['N_from']
    params['N_to'] = params['N_from'] + 100 if 'N_to' not in params.keys() else params['N_to']
    params['M_from'] = 0 if 'M_from' not in params.keys() else params['M_from']
    params['M_to'] = params['M_from'] + 10 if 'M_to' not in params.keys() else params['M_to']
    params['k_max'] = 5000 if 'k_max' not in params.keys() else params['k_max']
    params['combine'] = False if 'combine' not in params.keys() else params['combine']
    params['par'] = False if 'par' not in params.keys() else params['par']
    params['g_all'] = False if 'g_all' not in params.keys() else params['g_all']
    params['save_networks'] = False if 'save_networks' not in params.keys() else params['save_networks']
    params['logfile'] = 'log.txt' if 'logfile' not in params.keys() else params['logfile']
    return params

def default_missing_fields(fields):
    """Check fields dict. E.g. fields_dict = {'gender.demographic' : {'role': 'confounder', 'type': 'CATEGORY'}}. Further documentation on the fields_dict can be found at TODO.
    If missing, add default values role = 'confounder' and type = 'CATEGORY in any entry in fields dict.
    Parameters
    ----------
    fields : dict
        Dictionary populated with user input about fields, i.e. confunders and variables.
    Return
    ----------
    fields : dict
        Checked fields dict.
    """
    fields = {} if fields is None else fields
    for conf_dict in fields.values():
        conf_dict['role'] = 'confounder' if 'role' not in conf_dict.keys() else conf_dict['role']
        conf_dict['type'] = 'CATEGORY' if 'type' not in conf_dict.keys() else conf_dict['type']
    return fields

def verify_input(data, params, fields):
    data = default_missing_data(data)
    params = default_missing_params(params)
    fields = default_missing_fields(fields)
    if len(fields.keys()) == 0:
        assert params['g_all'] == True, 'If no fields are specified, g_all must be set to True. Otherwise, nothing to do here.'
    for key in fields.keys():
        print(key)
        fields[key]['type'] = Selectors.BlockType(fields[key]['type'])
        fields[key]['role'] = Selectors.Role(fields[key]['role'])
    assert params['g_all'] or len(list(data.keys())) > 0, 'If no confounders are specified, -g_all flag must be set to infer network from entire data instead.'
    assert params['N_from'] <= params['N_to'], 'N_from must be smaller than or equal to N_to'
    assert params['M_from'] <= params['M_to'], 'M_from must be smaller than or equal to M_to'
    for key in data.keys():
        if data[key]['tcga']:
            data[key] = {'ged': get_tcga_ged_name(key), 'pt': get_tcga_pt_name(key), 'sep': ',', 'tcga': True, 'tissue_type_field': data[key]['tissue_type_field'], 'tissue_type': data[key]['tissue_type']}
        else:
            assert data[key]['ged'] and data[key]['pt'], 'Specify ged (gene expression data) file name and pt (pheno type) file name for each cohort if tcga option is set to False.'

def get_tcga_ged_name(key):
    """Get string of TCGA gene expression file for cohort specified in key. TCGA files can be downloaded using download_tcga_cohorts.py at TODO.
    Parameters
    ----------
    key : str
        Study abbreviation of TCGA cohort.
    Return
    ----------
    : str
        File name of gene expression file corresponding to cohort specified in key.
    """
    assert key in [str(ident) for ident in list(Selectors.TCGACancerTypeSelector)]
    return f'TCGA-{key}.htseq_fpkm.tsv'

def get_tcga_pt_name(key):
    """Get string of TCGA pheno type file for cohort specified in key. TCGA files can be downloaded using download_tcga_cohorts.py at TODO.
    Parameters
    ----------
    key : str
        Study abbreviation of TCGA cohort.
    Return
    ----------
    : str
        File name of pheno type file corresponding to cohort specified in key.
    """
    assert key in [str(ident) for ident in list(Selectors.TCGACancerTypeSelector)]
    return f'TCGA-{key}.GDC_phenotype.tsv'
