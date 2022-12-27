import pandas as pd
import os
import Selectors

cohort = Selectors.CancerTypeSelector.BRCA

cwd = os.getcwd()
pam_path = os.path.join(cwd, 'data', 'brca_pam50')
pheno = Selectors.get_pheno_data(cohort)
pam = pd.read_csv(pam_path, sep='\t', header=0)
pam = pam.rename({'sample': 'submitter_id.samples'}, axis='columns')
pam = pam[pam.duplicated() == False]
pheno = pd.merge(pheno, pam, on="submitter_id.samples")
assert pheno[pheno.duplicated(subset='submitter_id.samples') == True].empty

expr = Selectors.get_expression_data(cohort)

print(set(pheno['PAM50']))
for subt in set(pheno['PAM50']):
    cur = pheno[pheno['PAM50'] == subt]
    cur.to_csv(os.path.join(cwd, 'data', 'TCGA-BRCA' + str(subt) + '.GDC_phenotype.tsv'), sep='\t')

    cur_e = expr[expr.index.isin(cur['submitter_id.samples'])]
    cur_e.to_csv(os.path.join(cwd, 'data', 'TCGA-BRCA'+str(subt)+'.htseq_fpkm.tsv'), sep='\t')