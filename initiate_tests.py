from test_suite import Selectors

data_dir = os.path.join(cwd, 'data')

print('start preprocessing...')
if not os.path.exists(data_dir):
    os.mkdir(data_dir)
    
# TODO: define cancer_types to be investigated. E.g.:
Selectors.download_TCGA_expression_data(Selectors.CancerTypeSelector.BLCA)
print('saved pd.DataFrame to csv file in ' + str(data_dir))

Selectors.download_known_tfs()
print('saved pd.DataFrame to csv file in ' + str(data_dir))
