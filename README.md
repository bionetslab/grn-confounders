# Investigating confounding in network inference - the confinspect package

Investigating confounding in gene regulatory networks (GRNs) and gene co-expression networks (GCNs) inferred from gene expression data.

The scripts and package provided in this repository allow users to investigate confounding in network inference using the methods ARACNe-AP, CEMiTool, GENIE3, GRNBoost2, and WGCNA, as presented in this [paper](TODO). The directory 'results_paper' contains all the results presented and jupyter-notebooks used in the paper. Further, programmers can import the confinspect package directly into their own scripts to implement customized algorithm wrappers and testing procedures.

Information on the output of the program can be found at the very bottom in 'Plot the results: output explained'.

## Prerequisites
Install [Python](https://www.python.org/downloads/) version 3.8, [R](https://www.r-project.org/) version 4.2, and [Java](https://www.oracle.com/java/technologies/downloads/) openjdk version "1.8.0_345". Install [Apache ANT](https://ant.apache.org/). Upgrade [pip](https://pypi.org/project/pip/).

Download this repository. In the directory results/plots/ you find the plots in our paper. In results/, you find the jupyter-notebook and all results used to create the plots.
#### Install confinspect as a local package from the source code
Move into the root directory of the repository, then run the following command:
```
pip install .
```
pip will install confinspect from the .toml file and the code in src/confinspect as a local package. Please note that the package might be removed upon disactivation of your environment, and re-installation via the command above might be necessary.
#### Install the required python packages
matplotlib==3.6.0
mpi4py==3.1.4
numpy==1.23.3
pandas==1.4.4
PyYAML==6.0
qnorm==0.8.1
scipy==1.9.1
seaborn==0.12.2
#### Download list of human known transcription factors and list of protein-coding genes
The list of protein-coding genes is required for all tests, the list of human known transcription factors is only required for testing of GRN inference methods. The list of protein-coding genes can be downloaded as follows:
```
python download_tcga_cohorts.py -pcgs
```
The list of human known transcription factors can be downloaded as follows:
```
python download_tcga_cohorts.py -tfs
```
List of protein-coding genes: http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/locus_groups/protein-coding_gene.txt

Seal RL, Braschi B, Gray K, Jones TEM, Tweedie S, Haim-Vilmovsky L, Bruford EA. Genenames.org: the HGNC resources in 2023. Nucleic Acids Res. PMID: 36243972 DOI: 10.1093/nar/gkac888

HGNC Database, HUGO Gene Nomenclature Committee (HGNC), European Molecular Biology Laboratory, European Bioinformatics Institute (EMBL-EBI), Wellcome Genome Campus, Hinxton, Cambridge CB10 1SD, United Kingdom www.genenames.org. Date of access: 01.10.2022

List of human known transcription factors: http://humantfs.ccbr.utoronto.ca/download/v_1.01/DatabaseExtract_v_1.01.csv

Lambert SA, Jolma A, Campitelli LF, Das PK, Yin Y, Albu M, Chen X, Taipale J, Hughes TR, Weirauch MT.(2018) The Human Transcription Factors. Cell. 172(4):650-665. doi: 10.1016/j.cell.2018.01.029. Review.

*IMPORTANT*

There are two options for users from here on: 1) use the scripts provided in the root directory of this repository to run tests with our predefined algorithm wrappers --> continue with the instructions right below. 2) Use the local package confinspect to implement your own algorithm wrapper (inherit from NetwrkInferenceWrapper.py) --> skip to section 'Option 2: Test your own methods'.
## Option 1: Run the tests using our predefined wrappers
#### Install the methods called by our predefined wrapper scripts
This step only has to be performed if you intend to use one of the predefined methods: ARACNe-AP, CEMiTool, GENIE3, GRNBoost2, or WGCNA. If you want to test other methods, please refer to the information given in section 'Option 2: Test your own methods'. 
##### ARACNe-AP
After downloading the [ARACNe-AP source code](https://github.com/califano-lab/ARACNe-AP), use the following command in the root directory of ARACNe-AP to build the ``jar`` executable:
```
ant main
```
The jar will be placed in ``dist/aracne.jar``. The documentation can be found in ``docs/index.html``.

For more information on using ARACNe-AP, check this [Readme](https://github.com/bionetslab/grn-confounders/tree/main/algorithms/ARACNe-AP)
##### GENIE3
Huynh-Thu V, Irrthum A, Wehenkel L, Geurts P (2010). “Inferring regulatory networks from expression data using tree-based methods.” PLoS ONE, 5(9), e12776. doi: 10.1371/journal.pone.0012776.

Aibar S, Bravo Gonzalez-Blas C, Moerman T, Huynh-Thu V, Imrichova H, Hulselmans G, Rambow F, Marine J, Geurts P, Aerts J, van den Oord J, Kalender Atak Z, Wouters J, Aerts S (2017). “SCENIC: Single-Cell Regulatory Network Inference And Clustering.” Nature Methods, 14, 1083-1086. doi: 10.1038/nmeth.4463.

Run the folloing commands in the R terminal:
```bash
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("GENIE3")
```   
 ##### WGCNA
Langfelder P, Horvath S (2008). “WGCNA: an R package for weighted correlation network analysis.” BMC Bioinformatics, 559. https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559.

Langfelder P, Horvath S (2012). “Fast R Functions for Robust Correlations and Hierarchical Clustering.” Journal of Statistical Software, 46(11), 1–17. https://www.jstatsoft.org/v46/i11/.

 Run the following commands in the R terminal:
```bash
 if (!require("WGCNA", quietly = TRUE))
   install.packages("WGCNA", repos = "http://cran.us.r-project.org")
```   
##### CEMiTool
Russo PdST, Ferreira GR, Cardozo LE, Burger MC, Arias-Carrasco R, Maruyama SR, Hirata TDC, Lima DS, Passos FM, Fukutani KF, Lever M, Silva JS, Maracaja-Coutinho V, Nakaya HTI (2018). “CEMiTool: a Bioconductor package for performing comprehensive modular co-expression analyses.” BMC Bioinformatics, 19(56), 1–13. doi: 10.1186/s12859-018-2053-1, https://doi.org/10.1186/s12859-018-2053-1.

Run the following commands in the R terminal:
```bash
 if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("CEMiTool")
```
##### GRNBoost2
Thomas Moerman, Sara Aibar Santos, Carmen Bravo González-Blas, Jaak Simm, Yves Moreau, Jan Aerts, Stein Aerts, GRNBoost2 and Arboreto: efficient and scalable inference of gene regulatory networks, Bioinformatics, Volume 35, Issue 12, June 2019, Pages 2159–2161, https://doi.org/10.1093/bioinformatics/bty916

Run the following command line to install the python arboreto package containing the GRNBoost2 method:
```
pip install arboreto
```
#### Test your own data
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
Make sure you are connected to the internet when running download_tcga_data.py. The script will download the data from [UCSC Xena](TODO), aggregate substages into superior stages and merge stages iii and iv, as described in our paper. Note that the columns 'gender.demographic', 'age_at_initial_pathologic_diagnosis', 'race.demographic', and 'tumor_stage.diagnoses' will not be renamed into 'sex', 'age', 'ethnicity', and 'stage', respectively, as in our tests. The pheno type file and gene expression file are saved to the data directory.

To download the LUSC and HNSC cohorts and the human known transcription factors, use the following command:
```
python download_tcga_cohorts.py -ct LUSC HNSC -tfs -pcgs
``` 
All files are saved to the data directory. The list of human known transcription factors is saved to 'regulators.csv'. The list of protein-coding genes is saved to 'protein-coding_gene.csv'.
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
python run_tests.py
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
##### Plot the results: output explained
*Jupyter-notebook used in our paper can be found in results/

Visualization of mean JIs obtained from comparing the networks inferred from the different blocks (see first part of the test protocol described in our paper): check whether the computed mean Jaccard Indices (JIs) are stored in results/JI. The path names resolve as follows:

JIs obtained from confounder-based partition _i_ of the cohort _cohort_ induced by the variable _variable_ inferred by the method _method_: cb\_{i}\_{method}\_{variable}\_{cohort}\_jaccInd.csv

JIs obtained from random partition _i_ of the cohort _cohort_ induced by the variable _variable_ inferred by the method _method_: rnd\_{i}\_{method}\_{variable}\_{cohort}\_jaccInd.csv

Modify the entries for method names, confounders, and cohorts in the jupyter-notebook visualize_mean_JIs_blocks.ipynb, then run the notebook.

Visualization of mean JIs obtained from comparing each network inferred from a single block with a network inferred from the whole data (see second part of the test protocol described in our paper): check whether the computed mean Jaccard Indices (JIs) are stored in results/JI. The path names resolve as follows:

JIs obtained from comparison of the network inferred from block _l_ belongin to the confounder-based partition _i_ of the cohort _cohort_ induced by the variable _variable_, with the network inferred from the entire data, both inferred by the method _method_: g\_all\_conf\_i\_{method}\_{confounder}\_{cohort}\_{l}\_jaccInd.csv

JIs obtained from comparison of the network inferred from block _l_ belongin to the random partition _i_ of the cohort _cohort_ induced by the variable _variable_, with the network inferred from the entire data, both inferred by the method _method_: g\_all\_rnd\_i\_{method}\_{confounder}\_{cohort}\_{l}\_jaccInd.csv

All lists of mean Jaccard Indices are comma-separated csv files with five columns: size intersection, size union, state (information about availability of edges in the compared networks), k, mean JI.

All networks are stored as comma-separated csv files with four columns: source gene, target gene, score, type. In case of undirected networks, the two connected genes are sorted lexicographically (source <= target).

log.txt is the default name for the log file. It contains information about block sizes, the chi^2 tests, and program run.

## License
GNU General Public License v3.0
