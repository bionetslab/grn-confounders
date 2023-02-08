# Investigating confounding in network inference

The scripts in this repository implement the test protocol described in TODO, using the confinspect package that was developed for this purpose. The directory 'results_paper' contains all results and jupyter-notebooks used in the paper. Detailed information on the confinspect package can be found at TODO.

### confinspect
confinspect is a Python package for investigating confounding by phenotypic variables in gene regulatory network (GRN) inference and gene co-expression network (GCN) inference. The core class of the package is the TestRunner, which allows users to organize tests involving different algorithms, data sets, and variables. The package further includes a collection of predefined wrappers for popular gene network inference algorithms (ARACNe-AP, CEMiTool, GENIE3, GRNBoost2, WGCNA), as well as tools for data preprocessing and selection. 

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
### Download this repository
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
##### Plot the results of the first part of the test
Check whether the computed mean Jaccard Indices are stored in results/JI, and the networks are stored in results/networks. Only the networks 
##### Plot the results of the second part of the test




To contribute to the development of confinspect, please refer to the CONTRIBUTING file for guidelines and instructions.
License


















## Installation


          

## Running the framework

After the tools are installed, clone the framework on your local machine.
Open terminal in the project folder and run the following command:

```python
python run_tests.py -ct BLCA -conf age -alg ARACNE -k 500 seq -n 100 -m 1

# -ct stands for disease (You may select: PCPG, GBM, COAD, BRCA, LUAD, PRAD, SKCM)
# -conf stands for confounders (You may select: sex, race, age, stage, type)
# -alg stands for algorithm (You may select: ARACNE, GENIE3, WGCNA, CEMI)
# -k stands for the top k edges of the network (You may select k value between 1 and 1000)
# -n stands for the number of partitions (You may select n between 1 and 1000)
# -m stands for random procedures. For GENIE3 select -m as 100. For other algorithms -m is 1.
```

NOTE: The user should be connected to the internet for downloading the gene dataset(TCGA). If the user is not connected to the network, the gene dataset has to be manually added by the user in the data folder. 

## Result Plots
                                           ARACNe result plot
![ARACNE(111)](https://user-images.githubusercontent.com/106863105/191177650-a46fb3c7-9194-4f4d-9622-052b8e4fcb8a.png)

                                            WGCNA result plot
![WGCNA_1](https://user-images.githubusercontent.com/106863105/191212093-54a66e60-4614-47cd-bc3e-52733aa6ec83.png)




## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.


## License
[MIT](https://choosealicense.com/licenses/mit/)
