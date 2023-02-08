# Investigating confounding in network inference

grn-confounders is a collection of scripts that allow users to investigate confounding in network inference. The scripts provided in this repository implement the test protocol presented in TODO. The directory 'results_paper' contains all results and jupyter-notebooks used in the paper. For more customized testing, programmers can import the confinspect package directly into their own scripts. Detailed information on the confinspect package can be found at TODO.

### Prerequisites
Install [Python](https://www.python.org/downloads/) version 3.8, [R](https://www.r-project.org/) version 4.2, and [Java](https://www.oracle.com/java/technologies/downloads/) openjdk version "1.8.0_345". Install [Apache ANT](https://ant.apache.org/). Upgrade pip.
##### Install python packages
confinspect==TODO
matplotlib==3.6.0
mpi4py==3.1.4
numpy==1.23.3
pandas==1.4.4
PyYAML==6.0
qnorm==0.8.1
scipy==1.9.1
seaborn==0.12.2
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
## Instructions on reproducing the tests presented in the paper
To reproduce the presented tests, follow the steps below:
##### Download this repository
After downloading this repository, check for the following directries and make sure the algorithms directory contains the caller script belonging to the predefined algrithm wrappers.

root/
  algorithms/
  config/
  data/
  results/
    JI/
    networks/
  temp/

##### Download TCGA data
Use the command line tool of the downloader script to download the gene expression data and pheno type data of a TCGA cohort, exemplarily the LUSC cohort. Make sure you are connected to the internet and have permission to access the internet on your device. The downloader will rename the fields ... and aggregate substages into the superior stages i, ii, and a joint categoriy iii and iv.
```
TODO downloader cmd
``` 
##### Run the tests
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
