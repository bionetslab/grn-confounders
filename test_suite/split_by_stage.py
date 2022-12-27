import pandas as pd
import os
import Selectors
cwd = os.getcwd()
#cohorts = [Selectors.CancerTypeSelector.KIRC, Selectors.CancerTypeSelector.KIRP, Selectors.CancerTypeSelector.LUSC, Selectors.CancerTypeSelector.BRCA]
cohorts = [Selectors.CancerTypeSelector.COAD]
for cohort in cohorts:
    #Selectors.get_expression_data(cohort)
    
    pheno = Selectors.get_pheno_data(cohort)
    pheno_field = 'tumor_stage.diagnoses'
    pheno = pheno[pheno[pheno_field] != 'stage x']
    pheno.loc[pheno['tumor_stage.diagnoses'].str.strip().isin(['stage ia', 'stage ib', 'stage ic']), pheno_field] = 'stage i'
    pheno.loc[pheno['tumor_stage.diagnoses'].str.strip().isin(['stage iia', 'stage iib', 'stage iic']), pheno_field] = 'stage ii'
    pheno.loc[pheno['tumor_stage.diagnoses'].str.strip().isin(['stage iiia', 'stage iiib', 'stage iiic', 'stage iv', 'stage iva', 'stage ivb', 'stage ivc']), pheno_field] = 'stage iii'

    for stage in set(pheno['tumor_stage.diagnoses']):
        cur = pheno[pheno['tumor_stage.diagnoses'] == stage]
        cur.to_csv(os.path.join(cwd, 'data', 'TCGA-' + str(cohort)  + str(stage) + '.GDC_phenotype.tsv'), sep='\t')

        #cur_e = expr[expr.index.isin(cur['submitter_id.samples'])]
        #cur_e.to_csv(os.path.join(cwd, 'data', 'TCGA-' + str(cohort) + str(stage)+'.htseq_fpkm.tsv'), sep='\t')
    
    
