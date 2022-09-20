# Investigating algorithmic bias in computational systems biology tools

This project is about assessing to what extent are high computational systems biology methods robust w.r.t data biases induced by the confounders such as age, sex, ethnicity.

## Systems Biology Tools

We focused on tools which uses gene data to create networks such as Gene Regulatory Network (Directed Network Graph)  and Gene Co-Expression Network (Undirected Network Graph). 

- **Gene Regulatory Network**

  - GENIE3 [repository](https://bioconductor.org/packages/release/bioc/html/GENIE3.html)
  - ARACNe [repository](https://bio.tools/aracne)

- **Gene Co-expression Network**
  - WGCNA [repository](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/)
  - CEMiTool [repository](https://www.bioconductor.org/packages/release/bioc/html/CEMiTool.html)

## Installation

***Prerequisites***: Install lastest versions of [Python](https://www.python.org/downloads/), [R](https://www.r-project.org/) and [Java](https://www.oracle.com/java/technologies/downloads/) based on your OS.

**Installation of GRN Tools**

   GENIE3 Tool: Run these commands in the R terminal 
```bash
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("GENIE3")
```    
ARACNe-AP Tool:  ``ARACNe-AP`` requires JDK > 1.8 and ANT. Use the following command in the repository root directory to build the ``jar`` and documentation:

```
ant main
```

The jar will be placed in ``dist/aracne.jar``. The documentation can be found in ``docs/index.html``.

For more information on using ARACNe-AP check this [Readme](https://github.com/bionetslab/grn-confounders/tree/main/algorithms/ARACNe-AP) 

**Installation of Co-Expression Tools**

 WGCNA Tool: Run these commands in the R terminal 
```bash
 if (!require("WGCNA", quietly = TRUE))
   install.packages("WGCNA", repos = "http://cran.us.r-project.org")
```    
 CEMi Tool: Run these commands in the R terminal 
```bash
 if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("CEMiTool")
```
          

## Running the framework

After the tools are installed, clone the framework on your local machine.
Open terminal in the project folder and run the command:

```python
python run_tests.py -ct BLCA -conf age -alg ARACNE -k 500 seq -n 100 -m 2 

# -ct stands for disease (You may select: PCPG, GBM, COAD, BRCA, LUAD, PRAD, SKCM)
# -conf stands for confounders (You may select: SEX, RACE, AGE, STAGE, TYPE)
# -alg stands for algorithm (You may select: ARACNE, GENIE3, WGCNA, CEMI)
# -k stands for the top k edges of the network (You may select k value between 1 and 1000)

```

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.


## License
[MIT](https://choosealicense.com/licenses/mit/)
