import Selectors
from test_runner import TestRunner
import pandas as pd
import numpy as np
import os

def test_conf_partition(conf_partition):
    assert all([len(el) == len(set(el)) for el in conf_partition])
    conf_partition = [set(el) for el in conf_partition]
    assert conf_partition.pop().intersection(*conf_partition) == set()

def test_rnd_partitions(rnd):
    # blocks do not overlap
    rnd = [[set(block) for block in cur] for cur in rnd]
    for part in rnd:
        assert part.pop().intersection(*part) == set()

if __name__ == '__main__':
    test_runner = TestRunner({'COAD':{'ged':'TCGA-COAD.htseq_fpkm.tsv', 'pt':'TCGA-COAD.GDC_phenotype.tsv'}}, 
            {'age_at_initial_pathologic_diagnosis':Selectors.BlockType.QUARTILE, 'tumor_stage.diagnoses':Selectors.BlockType.CATEGORY}, ['WGCNA'], 0, 1, 0, 1, 100)

    for ct_sel in test_runner.cancer_type_selectors:
        for conf_sel in test_runner.confounder_selectors:
            test_conf_partition(list(test_runner.conf_partitions[ct_sel][conf_sel].values()))
            test_rnd_partitions(test_runner.rnd_partitions)

    #res = pd.concat([pd.DataFrame({'mean JI':expr, 'k':range(10, 1000, 100)}) for expr in [range(10, 1000, 100), range(10, 1000, 100)]])
    #print(res)
    test_runner.run_on_cancer_types_confounders()