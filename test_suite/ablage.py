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

""" unter visualize JI, gleich eingerückt:
    for bl in block_ids:
        visualize...
        sign_conf_rnd = 0
        for k in rep_k:
            rnd = networks_blocks_rnd[bl][networks_blocks_rnd[bl]['k'] == k]
            conf = networks_blocks_conf[bl][networks_blocks_conf[bl]['k'] == k]
            wilcox_ = wilcoxon(conf['mean JI'], rnd['mean JI'], alternative='less', correction=True)
            if wilcox_.pvalue < 0.05:
                sign_conf_rnd += 1
        print(['One-sided wilcoxon test on ' + str(bl) + ' network and corresponding random network for ' +
                    str(sign_conf_rnd/len(rep_k)) + ' of all tested k.'])

    for bl0, bl1 in itt.product(block_ids, 2):
        net0 = networks_blocks_conf[bl0]
        net1 = networks_blocks_conf[bl1]
        sign_conf_conf = 0
        for k in rep_k:
            block0_conf = net0[net0['k'] == k]
            block1_conf = net1[net1['k'] == k]
            try:
                wilcox = wilcoxon(block0_conf['mean JI'],block1_conf['mean JI'], alternative='less', correction=True)
                if wilcox.pvalue < 0.05:
                    sign_conf_conf += 1
            except:
                assert all(x == 0.0 for x in (conf_0['mean JI']-conf_1['mean JI']))
                continue
        
        wilcox_fractions = pd.concat([wilcox_fractions, pd.DataFrame({'cohort':ct_sel, 'confounder/variable':str(conf_sel), 
                        'method':str(alg_sel), 'bl0':str(bl0), 'bl1':str(bl1), 'frac':sign_conf_conf/len(rep_k)})])
    return wilcox_fractions

In erste bzw. zweite Funktion:
    wilcox_fractions = pd.DataFrame(columns=['cohort', 'confounder/variable','method', 'bl0', 'bl1', 'frac'])

    wilcox_fractions = pd.DataFrame(columns=['cohort', 'confounder/variable','method', 'bl0', 'bl1', 'frac'])

"""