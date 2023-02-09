# Investigating confounding in network inference - the confinspect package

Investigating confounding in gene regulatory networks (GRNs) and gene co-expression networks (GCNs) inferred from gene expression data.

The scripts and package provided in this repository allow users to investigate confounding in network inference, as presented in this [paper](TODO). The directory 'results_paper' contains all the presented results and jupyter-notebooks used in the paper. Further, programmers can import the confinspect package directly into their own scripts for more customized testing.

There are two options for users: 1) use the scripts provided in the root directory of this repository to run tests with our predefined algorithm wrappers --> continue with the instructions right below. 2) Use the local package confinspect to implement your own algorithm wrapper (inherit from NetwrkInferenceWrapper.py) --> skip to section TODO.

## Prerequisites
Install [Python](https://www.python.org/downloads/) version 3.8, [R](https://www.r-project.org/) version 4.2, and [Java](https://www.oracle.com/java/technologies/downloads/) openjdk version "1.8.0_345". Install [Apache ANT](https://ant.apache.org/). Upgrade pip.

Download this repository. The directory result_paper can be omitted, unles you wish to visualize the results from our paper.
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
### 1) Run the tests using our predefined wrappers
To use the predefined wrappers, i.e. using one of the methods ARACNe-AP, CEMiTool, GENIE3, GRNBoost2, or WGCNA, you can use the runner script run_tests.py. The runner script takes all input parameters from the config files provided in the config directory. The config files are organized as follows:
    data.yml: information abut the gene expression data and pheno type data to be used. Example:
    HNSC: # arbitrary string identifier of the cohort
       tcga: !!bool False # iff the data was downloaded using the download_tcga_data.py TODO. Default is False
       ged: 'ged_file_name.csv' # name of the gene expression file, located in the data directory
       pt: 'pt_file_name.csv' # name of the pheno type file, located in the data directory
       sep: ',' # separator used in both the files at ged and pt, default is ','
       tissue_type_field: sample_type.samples # name of the column (field) in the pt file to be used to filter the samples. Must exist in the pt file.
       tissue_type: Primary Tumor # keep only samples that have 'Primary Tumor' in the 'sample_type.samples' column. Must exist in the pt file.
    LUSC:
        ...and so on
        
    fields.yml: information about the confounders
    gender.demographic: # name of the column in the pt file to be tested for confounding (i.e. used to induce a partition)
       role: CONFOUNDER # CONFOUNDER | VARIABLE - perform a chi2 test for dependency of all variables with all confounders
       type: CATEGORY # CATEGORY | QUARTILE - how to induce the partition. CATEGORY: split samples into blocks based on the exact expression in the column used for partitioning (e.g. male, female). QUARTILE: form a block of each of the upper and lower quartile accoridng to the values in the column used for partitioning. Column values must be numerical (e.g. upper and lower quartile of a numerical age column)
    tumor_stage.diagnoses:
       role: VARIABLE
       type: CATEGORY
       
    params.yml: all other parameters
    algorithms: # list methods to be tested below, in capital letters.
    - CEMITOOL
    - WGCNA
    N_from: 0 # starting index of random partitions per (cohort, confounder)
    N_to: 10 # stopping index of random partitions per (cohort, confounder)
    M_from: 0 # starting index 
    M_to: 10
    k: 5000
    combine: !!bool False
    par: !!bool False
    g_all: !!bool True
    save_networks: !!bool False
    logfile: log.txt
    
   
To reproduce the presented tests, you can use the scripts download_tcga_data.py and run_tests.py. The actual runner script, run_tests.py, is instructed by the config files data.yml, fields.yml, and params.yml in the config directory. You can use the files as-is, or modify them to run tests on your own data and confounders. Follow the steps below:
##### Download TCGA data
Use the command line tool of the downloader script to download the gene expression data and pheno type data of a TCGA cohort, exemplarily the LUSC cohort. Make sure you are connected to the internet and have permission to access the internet on your device. The downloader will rename the fields ... and aggregate substages into the superior stages i, ii, and a joint categoriy iii and iv. The data is stored in the data directory.
```
TODO downloader cmd
``` 
Use the ... script to start the actual tests. The config files provided in ... will instruct the script to test for confounding of the ... confounders in the LUSC cohort using the methods ... Further, test for an effect of the stage variable, and compute a $\chi^{2}$ test of the stage variable. 
```
TODO run_tests cmd
``` 
##### Plot the results
Check whether the computed mean Jaccard Indices are stored in results/JI. The path names resolve as follows: TODO

Run the jupyter-notebook ... located in the root directory to create line plots and swarm plots of the results.

Then check whether the mean Jaccard Indices of the cmparison with the networks inferred from the entire data are stored in results/JI. The path names resolve as follows: TODO

Run the jupyter-notebook ... located in the root directory to create line plots of the results.

## Instructions on running tests on your own data
TODO

## Instructions on running tests with new methods
To define custom algorithm wrappers, the user is referred to the script 'dummy_custom_' TODO, implementing the stub of a new method wrapper inheriting from NetwrokInferenceWrapper.py. For further instructions, please visit confinspect TODO and read the documentation and README.


## License
GNU GPL
