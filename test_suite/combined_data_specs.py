import pandas as pd
import os
import Selectors
"""
types = [Selectors.CancerTypeSelector.STAD, Selectors.CancerTypeSelector.COAD, Selectors.CancerTypeSelector.READ]

expr = {sel: Selectors.get_expression_data(sel) for sel in types}
pheno = {sel: Selectors.get_pheno_data(sel) for sel in types}
ct_sel = '-'.join([str(el) for el in types])
expr = pd.concat(expr.values())
pheno = pd.concat(pheno.values())

out_path_expr = f'{str(ct_sel)}_combined_data_specs.txt'
with open(out_path_expr, 'w') as k:
    k.write('COHORT: ' + str(ct_sel)+'\n')

    k.write('pheno data before aligning expr and pheno: ' + str(pheno.shape)+'\n')
    k.write('Remove genes with std=0 and align expression data and phenotype data on samples for cohort ' + str(ct_sel) + '...\n')
    keep = pheno['submitter_id.samples'].isin(expr.index)
    pheno = pheno[keep]
    pheno = pheno[pheno['submitter_id.samples'].isin(expr.index)]
    samples = pheno['submitter_id.samples']            
    expr = expr.loc[samples]
    #k.write('Remove gene where standard deviation of expression data is 0 for cohort ' + str(ct_sel) + '...\n')
    expr = expr.loc[:, (expr.std() != 0)]
    pheno = pheno[pheno['submitter_id.samples'].isin(expr.index)]
    k.write('---\n')
    k.write('final shapes: \n')
    k.write('expr shape: ' + str(expr.shape) +'\n')
    k.write('pheno shape: ' + str(pheno.shape)+'\n')
    #pheno.to_csv('combined_HNSC_LUSC.tsv', sep='\t')
    for conf_sel in list(Selectors.ConfounderSelector):
        k.write('CONFOUNDER: ' + str(conf_sel)+'\n')
        conf_partition, printer = Selectors.get_conf_partition(pheno, conf_sel, 0)
        i = 0
        for el in printer:
            #k.write(f'{str(blocks[i])}: ' + str(len(block))+'\n')
            k.write(str(el[1]) + ': '+str(len(el[0])) + '\n')
            for cohort in types:
                cur = pheno[pheno['submitter_id.samples'].isin(el[0])]
                cur = cur[cur['cohort'] == str(cohort)]
                k.write('of which '+ str(len(cur)) + ' are from cohort ' + str(cohort) + '\n')
            i += 1

    k.write('-----------------------------------\n')

"""


expr = Selectors.get_expression_data(Selectors.CancerTypeSelector.LUSC)
pheno = Selectors.get_pheno_data(Selectors.CancerTypeSelector.LUSC)

print('Align expression data and phenotype data on samples for cohort ' )
keep = pheno['submitter_id.samples'].isin(expr.index)
pheno = pheno[keep]
pheno = pheno[pheno['submitter_id.samples'].isin(expr.index)]
samples = pheno['submitter_id.samples']            
expr = expr.loc[samples]
print('Remove genes where standard deviation of expression data is 0 for cohort ')
expr = expr.loc[:, (expr.std() != 0)]
pheno = pheno[pheno['submitter_id.samples'].isin(expr.index)]

conf = Selectors.get_conf_partition(pheno, Selectors.ConfounderSelector.RACE, 0)
rnd = Selectors.get_n_random_partitions(0, 100, pheno['submitter_id.samples'], conf, Selectors.CancerTypeSelector.LUSC, Selectors.ConfounderSelector.RACE)

for block in conf:
    print(len(block))


for part in rnd:
    for block in part:
        print(len(block))
    b_1 = part[0]
    b_2 = part[1]
    #print(set(b_1).intersection(set(b_2)))






