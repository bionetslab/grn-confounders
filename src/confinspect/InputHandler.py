import argparse
import os
from confinspect import Selectors

def setup_directories():
    """Set up directory structure needed by the TestRunner and NetworkInferenceWrapper. Root is os.getcwd().
    """
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
    assert len(data.keys()) > 0, 'specify at least one cohort'
    for ct_dict in data.values():
        ct_dict['tcga'] = False if 'tcga' not in ct_dict.keys() else ct_dict['tcga']
        if 'sep' not in ct_dict.keys():
            print('Setting sep separator to \',\'...')
            ct_dict['sep'] = '\t'
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
    params = {} if params is None else params
    if 'algorithms' not in params.keys():
        print('No algorithm specified in params[\'algorithms\']. Please add custom algorithm after instantiating TestRunner by calling add_custom_algorithm before calling induce_partitions.')
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
    assert len(fields.keys()) > 0, 'specify at least one field'
    for conf_dict in fields.values():
        conf_dict['role'] = 'confounder' if 'role' not in conf_dict.keys() else conf_dict['role']
        conf_dict['type'] = 'CATEGORY' if 'type' not in conf_dict.keys() else conf_dict['type']
    return fields

def dump_config(data, params, fields, logger=None):
    """Dump defaulted config files to inform the user about the final configurations of his test runs.
    Parameters
    ----------
    args: argparse.Namespace
        Namespace object populated with user command line arguments.
    """
    cwd = os.getcwd()
    stamp = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    if not os.path.exists(config_path):
        os.mkdir(config_path)
    with open(os.path.join(cwd, 'config', f'data_def_{stamp}.yml'), 'w') as f:
        yaml.dump(data, f)
    with open(os.path.join(cwd, 'config', f'fields_def_{stamp}.yml'), 'w') as f:
        yaml.dump(fields, f)
    with open(os.path.join(cwd, 'config', f'params_def_{stamp}.yml'), 'w') as f:
        yaml.dump(params, f)
    if logger:
        logger.info(f'Defaulted incomplete config files. Modified files were dumped to data_def_{stamp}.yml, params_def_{stamp}.yml,\
             fields_def_{stamp}.yml\ in config directory.')

def verify_input(data_, params_, fields_, logger=None):
    """Check and default data, params, and fields dictionaries. The dictionaries ill contain all necessary entries afterwards.
    Parameters
    ----------
    data : dict
        Dictionary with data identifiers. Keys and types: {cohort_name: {'tcga': bool, 'ged': str, 'pt': str, 'sep': str, 'tissue_type_field': str, 'tissue_type': str}}
    params : dict
        Dictionary with parameter identifiers. Keys and types: {'algorithms': list, 'N_from': int, 'N_to': int, 'M_from': int, 'M_to': int, 'k_max': int, 'combine': bool, 'par': bool, 'g_all': bool, 'save_networks': bool, 'logfile': str}
    fields : dict
        Dictionary with parameter identifiers. Keys and types: {variable_name : {'role': 'confounder' | 'variable', 'type': 'CATEGORY' | 'QUARTILE'}}
    Return
    ----------
    checked and defaulted data, params, fields dictionaries
    """
    data = default_missing_data(data_)
    params = default_missing_params(params_)
    fields = default_missing_fields(fields_)
    if len(fields.keys()) == 0:
        assert params['g_all'] == True, 'If no fields are specified, g_all must be set to True. Otherwise, nothing to do here.'
    assert params['g_all'] or len(list(data.keys())) > 0, 'If no confounders are specified, -g_all flag must be set to infer network from entire data instead.'
    assert params['N_from'] <= params['N_to'], 'N_from must be smaller than or equal to N_to'
    assert params['M_from'] <= params['M_to'], 'M_from must be smaller than or equal to M_to'
    for key in data.keys():
        if data[key]['tcga']:
            data[key] = {'ged': get_tcga_ged_name(key), 'pt': get_tcga_pt_name(key), 'sep': '\t', 'tcga': True, 'tissue_type_field': data[key]['tissue_type_field'], 'tissue_type': data[key]['tissue_type']}
        else:
            assert data[key]['ged'] and data[key]['pt'], 'Specify ged (gene expression data) file name and pt (pheno type) file name for each cohort if tcga option is set to False.'
    if not (data == data_) and (fields == fields_) and (params == params_):
        dump_config(data, params, fields, logger)

    # change string identifiers into Selectors
    for key in fields.keys():
        assert str(fields[key]['type']) in [str(sel) for sel in list(Selectors.BlockType) if sel != Selectors.BlockType.ALL], f'type of {key} confounder must be one of ' + ','.join([str(sel) for sel in list(Selectors.BlockType) if sel != Selectors.BlockType.ALL])
        assert str(fields[key]['role']) in [str(sel) for sel in list(Selectors.Role)], f'role of {key} confounder must be one of ' + ','.join([str(sel) for sel in list(Selectors.Role)])
        fields[key]['type'] = Selectors.BlockType(fields[key]['type'])
        fields[key]['role'] = Selectors.Role(fields[key]['role'])
    return data, params, fields

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
