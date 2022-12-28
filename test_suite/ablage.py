"""
def _new_str(self):
    return self.value

def dynamicConfounderSelector(values):
    
    dynamic_enum = Enum('ConfounderSelector', values)
    dynamic_enum.__str__ = funcType(_new_str, dynamic_enum, dynamicConfounderSelector)
    return dynamic_enum

def dynamicCancerTypeSelector(values):
    #Dynamic enum specifying the cohorts to be examined. Is initialized based on the user input.
    #Parameters
    #----------
    #values : dict
        #Dictionary of enum elements with string values.

    #Returns
    #----------
    #dynamic_enum : Enum
        #Custom enum with elements that represent the cohorts to be examined.
    
    dynamic_enum = Enum('CancerTypeSelector', values)
    dynamic_enum.__str__ = funcType(_new_str, dynamic_enum, dynamicCancerTypeSelector)
    return dynamic_enum

def parse_input_file():
    args = {}
    args['combine'] = False
    f = open(os.path.join(os.getcwd(), "input_parameters.txt"))
    input_data = f.readlines()
    for line in input_data:
        key, value = line.split(":")
        if key == 'combine'.strip():
            args[key.strip()] = True
        elif key.strip() in ['N_from', 'N_to', 'M_from', 'M_to', 'k']:
            args[key.strip()] = int(value.strip())
        elif key == 'mode'.strip():
            assert value.strip() == 'par' or value.strip() == 'seq'
            args[key.strip()] = value.strip()
        else:
            args[key.strip()] = value.strip()
    f.close()
    assert all([el in args.keys() for el in ['ct', 'conf', 'alg', 'N_from', 'N_to', 'M_from', 'M_to', 'k', 'mode']])
    args = access_with_dot(args)
    return args

class access_with_dot(dict):
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

"""

