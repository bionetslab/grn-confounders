import pandas as pd
import os
import Selectors

#os.chdir(os.path.join('BIONETS_code'))
for ct_sel in list(Selectors.CancerTypeSelector):
    pheno = Selectors.get_pheno_data(ct_sel)
    #out_path_pheno = f'{str(ct_sel)}_pheno_specs.txt'
    #with open(out_path_pheno, 'w') as f:
    #    f.write('COHORT: ' + str(ct_sel)+'\n')
    #    f.write('Pheno data shape: ' + str(pheno.shape)+'\n')
    #    for conf_sel in list(Selectors.ConfounderSelector):
    #        f.write('CONFOUNDER: ' + str(conf_sel)+'\n')
    #        conf_partition, blocks = Selectors.get_conf_partition(pheno, conf_sel)
    #        i = 0
    #        for block in conf_partition:
    #            f.write(f'{str(blocks[i]}: ' + str(len(block))+'\n')
    #            i+=1

            #f.write('Filter Primary Tumor samples in pheno data for cohort ' + str(ct_sel) + '...\n')
            #pheno =  pheno[pheno['sample_type.samples'] == 'Primary Tumor']
            #conf_partition = Selectors.get_conf_partition(pheno, conf_sel)
            #i = 0
            #for block in conf_partition:
            #    f.write(f'Block nb {i}: ' + str(len(block))+'\n')
            #    i += 1

        #f.write('-------------------\n')

    expr = Selectors.get_expression_data(ct_sel)
    out_path_expr = f'{str(ct_sel)}_data_specs.txt'
    with open(out_path_expr, 'w') as k:
        k.write('COHORT: ' + str(ct_sel)+'\n')
    #    k.write('expr data shape before prp: ' + str(expr.shape)+'\n')
    #    k.write('---\n')

    #    k.write('Remove version identifiers from gene symbols in expression data for cohort ' + str(ct_sel) + '...\n')
    #    expr.columns = expr.columns.str.split('.').str[0].tolist()
    #    
    #    k.write('Only leave protein-coding genes in expression data set for cohort ' + str(ct_sel) + '...\n')
    #    pcg = pd.read_csv(os.path.join('protein-coding_gene.csv'))
    #    k.write('p-c genes shape: ' + str(pcg.shape)+'\n')
    #    mask = expr.columns.intersection(pcg['ensembl_gene_id'].values)
    #    expr = expr[mask]
    #    k.write('expr data shape after pcg-filtering: ' + str(expr.shape)+'\n')
    #    k.write('---\n')

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
        for conf_sel in list(Selectors.ConfounderSelector):
            k.write('CONFOUNDER: ' + str(conf_sel)+'\n')
            conf_partition, blocks = Selectors.get_conf_partition(pheno, conf_sel)
            i = 0
            for block in conf_partition:
                k.write(f'{str(blocks[i])}: ' + str(len(block))+'\n')
                i += 1

        k.write('-----------------------------------\n')


    



    
