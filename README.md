# Investigating confounding in network inference - the confinspect package

Investigating confounding in gene regulatory networks (GRNs) and gene co-expression networks (GCNs) inferred from gene expression data.

The scripts and package provided in this repository allow users to investigate confounding in network inference, as presented in this [paper](TODO). The directory 'results_paper' contains all the presented results and jupyter-notebooks used in the paper. Further, programmers can import the confinspect package directly into their own scripts for more customized testing.

## Prerequisites
Install [Python](https://www.python.org/downloads/) version 3.8, [R](https://www.r-project.org/) version 4.2, and [Java](https://www.oracle.com/java/technologies/downloads/) openjdk version "1.8.0_345". Install [Apache ANT](https://ant.apache.org/). Upgrade [pip](https://pypi.org/project/pip/).

Download this repository. The directory results_paper can be omitted, unless you only wish to visualize the results from our paper.
#### Install confinspect as a local package from the source code src/confinspect/. Move into the root directory of the repository, then run the following command:
```
pip install .
```
#### Install the required python packages
matplotlib==3.6.0
mpi4py==3.1.4
numpy==1.23.3
pandas==1.4.4
PyYAML==6.0
qnorm==0.8.1
scipy==1.9.1
seaborn==0.12.2

*IMPORTANT*

There are two options for users from here on: 1) use the scripts provided in the root directory of this repository to run tests with our predefined algorithm wrappers --> continue with the instructions right below. 2) Use the local package confinspect to implement your own algorithm wrapper (inherit from NetwrkInferenceWrapper.py) --> skip to section 'Option 2: Test your own methods'.
## Option 1: Run the tests using our predefined wrappers
#### Install the methods called by our predefined wrapper scripts
This step only has to be performed if you intend to use one of the predefined methods: ARACNe-AP, CEMiTool, GENIE3, GRNBoost2, or WGCNA. If you want to test other methods, please refer to the information given in section TODO. 
##### ARACNe-AP
After downloading the [ARACNe-AP source code](https://github.com/califano-lab/ARACNe-AP), use the following command in the root directory of ARACNe-AP to build the ``jar`` executable:
```
ant main
```
The jar will be placed in ``dist/aracne.jar``. The documentation can be found in ``docs/index.html``.

For more information on using ARACNe-AP, check this [Readme](https://github.com/bionetslab/grn-confounders/tree/main/algorithms/ARACNe-AP)
##### GENIE3
Run the folloing commands in the R terminal:
```bash
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("GENIE3")
```   
 ##### WGCNA
 Run the following commands in the R terminal:
```bash
 if (!require("WGCNA", quietly = TRUE))
   install.packages("WGCNA", repos = "http://cran.us.r-project.org")
```   
##### CEMiTool
Run the following commands in the R terminal:
```bash
 if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("CEMiTool")
```
##### GRNBoost2
Run the following command line to install the python arboreto package containing the GRNBoost2 method:
```
pip install arboreto
```
#### Run tests on your own data
(If you want to reproduce the tests from our paper, please go to the next section)

To use the predefined wrappers, i.e. using one of the methods ARACNe-AP, CEMiTool, GENIE3, GRNBoost2, or WGCNA, you can use the runner script run_tests.py. The runner script takes all input parameters from the config files provided in the config directory. The config files are organized as follows:

    data.yml: information abut the gene expression data and pheno type data to be used. Example:
    HNSC: # arbitrary string identifier of the cohort
       ged: 'ged_file_name.csv' # name of the gene expression file, located in the data directory
       pt: 'pt_file_name.csv' # name of the pheno type file, located in the data directory
       sep: ',' # separator used in both the files at ged and pt, default is ','
       tissue_type_field: sample_type.samples # name of the column (field) in the pt file to be used to filter the samples. Must exist in the pt file.
       tissue_type: Primary Tumor # keep only samples that have 'Primary Tumor' in the 'sample_type.samples' column. Must exist in the pt file.
    LUSC:
        ...and so on
        
    fields.yml: information about the confounders. We distinguish between variables and confounders ('special variables'). The chi2-tests are performed for each variable with all confounders. The results of the chi2 tests are printed to the log file.
    gender.demographic: # name of the column in the pt file to be tested for confounding (i.e. used to induce a partition)
       role: CONFOUNDER # CONFOUNDER | VARIABLE - perform a chi2 test for dependency for each variable with all confounders
       type: CATEGORY # CATEGORY | QUARTILE - how to induce the partition. CATEGORY: split samples into blocks based on the exact expression in the column used for partitioning (e.g. male, female). QUARTILE: form a block of each of the upper and lower quartile accoridng to the values in the column used for partitioning. Column values must be numerical (e.g. upper and lower quartile of a numerical age column)
    tumor_stage.diagnoses:
       ... and so on
       
    params.yml: all other parameters
    algorithms: # list methods to be tested below, in capital letters.
    - CEMITOOL
    - WGCNA
    N_from: 0 # starting index of random partitions per (cohort, confounder), default is 0
    N_to: 10 # stopping index of random partitions per (cohort, confounder), default is 100
    M_from: 0 # starting index of the copies of the confounder-induced partition per (cohort, confounder), default is 0
    M_to: 10 # stopping index of the copies of the confounder-induced partition per (cohort, confounder), default is 10
    k: 5000 # Compute mean Jaccard  Index of top k edges, k in range(10, k,100), default is 5000
    combine: !!bool False # combine all cohorts specified in data.yml into one cohort and test for the effect of 'cohort name' as confounder, default is False
    par: !!bool False # use MPI parallelization, default is False
    g_all: !!bool False # infer networks from entire data set, for overlapping indices of range(N_from, N_to) and range(M_from, M_to), default is False
    save_networks: !!bool False # save networks from which the mean Jaccard Indices were computed, default is False (high storage consumption)
    logfile: log.txt # print info log to this file (in 'w', i.e. overwrite mode), default is log.txt
    
You can modify the config files such that you can run your own tests, or use the config files given in the section below, to...
#### Reproduce the tests from the paper
To reproduce the presented tests, you can use the scripts download_tcga_data.py before running run_tests.py. All scripts will use the TCGA study abbreviations we use in our paper to address the different cohorts.
##### Download TCGA data
Make sure you are connected to the internet when running download_tcga_data.py. The script will download the data from [UCSC Xena](TODO), aggregate substages into superior stages and merge stages iii and iv, as described in our paper. Further, the columns 'gender.demographic', 'age_at_initial_pathologic_diagnosis', 'race.demographic', and 'tumor_stage.diagnoses' will be renamed into 'sex', 'age', 'ethnicity', and 'stage', respectively, as in our tests. The pheno type file and gene expression file are saved to the data directory.

To download the LUSC and HNSC cohorts, for example, use the following command:
```
TODO downloader cmd
``` 
##### Set up the config files

    data.yml: 
    HNSC:
       tcga: !!bool True
       tissue_type_field: sample_type.samples 
       tissue_type: Primary Tumor
    LUSC:
       tcga: !!bool True
       tissue_type_field: sample_type.samples 
       tissue_type: Primary Tumor
       
    fields.yml: 
    age: # name of the column in the pt file to be tested for confounding (i.e. used to induce a partition)
       role: CONFOUNDER
       type: QUARTILE
    sex:
       role: CONFOUNDER
       type: CATEGORY
    ethnicity:
       role: CONFOUNDER
       type: CATEGORY
    stage:
       role: VARIABLE
       type: CATEGORY       
       
    params.yml:
    algorithms:
    - ARACNE
    - CEMITOOL
    - GENIE3
    - GRNBOOST2
    - WGCNA
    N_from: 0 
    N_to: 100 
    M_from: 0 
    M_to: 10 
    k: 5000 
    par: !!bool True 
    g_all: !!bool True

These config files can be found in the config directory with the appendix '\_tcga'. Remove the appendix or copy the lines into the data.yml, fields.yml, params.yml files.
##### Start the tests
Use the run_tests.py script to start the actual tests. The config files from above will instruct the script.
```
TODO run_tests cmd
``` 
## Option 2: Test your own methods
To define custom algorithm wrappers, the user can use the following code stub as a blue print:
```python
import confinspect
from confinspect import NetworkInferenceWrapper, TestRunner, Selectors
import InputHandler
import os

class CustomWrapper(NetworkInferenceWrapper.NetworkInferenceWrapper):

    def _infer_network(self, expression_data, rank):
        print('infer network from a single block!')
        pass

    def _get_top_k_edges(self, i, k):
        print('Get top k edges, with each edge being a tuple of nodes, and don\'t forget to order the nodes alphabetically, if the \
            edges are undirected!')
        pass

if __name__ == "__main__":
    InputHandler.setup_directories()
    wrap = CustomWrapper()
    data = {'HNSC': {'tcga': True, 'ged': f'TCGA-HNSC.htseq_fpkm.tsv', 'pt': f'TCGA-HNSC.GDC_phenotype.tsv', 'sep': ',', 'tissue_type_field': None, 'tissue_type': None}}
    fields = {'gender.demographic': {'role': Selectors.Role.CONFOUNDER, 'type': Selectors.BlockType.CATEGORY}}
    params = {'N_from': 0, 'N_to': 1, 'M_from': 0, 'M_to': 1, 'k_max': 100, 'save_networks': False, 'g_all': True, 'logfile': 'log.txt'}
    InputHandler.verify_input(data, params, fields)
    tr = TestRunner.TestRunner(data, fields, params)
    tr.add_custom_algorithm(wrap, 'name')
    tr.induce_partitions()

    print('This will throw an error unless the CustomWrapper gets implemented properly.')
    tr.run_all()
```
##### Plot the results
Check whether the computed mean Jaccard Indices are stored in results/JI. The path names resolve as follows: TODO

Run the jupyter-notebook ... located in the root directory to create line plots and swarm plots of the results.

Then check whether the mean Jaccard Indices of the cmparison with the networks inferred from the entire data are stored in results/JI. The path names resolve as follows: TODO

Run the jupyter-notebook ... located in the root directory to create line plots of the results.

## License
GNU GPL
