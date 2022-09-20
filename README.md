# Investigating algorithmic bias in computational systems biology tools

This project is about assessing to what extent are high computational systems biology methods robust w.r.t data biases induced by the confounders such as age, sex, ethnicity. Moreover, investigating the effects the confounders have on the Gene Regulatory Network and Co-expression network. Which results in the conclusion that, if there are major differences induced by a confounder, the data for these groups should be studied or at looked at separately.

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

**Installation of Co-Expression Network Tools**

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
![Screenshot (111)](https://user-images.githubusercontent.com/106863105/191177650-a46fb3c7-9194-4f4d-9622-052b8e4fcb8a.png)

                                            WGCNA result plot
![WGCNA_1](https://user-images.githubusercontent.com/106863105/191177823-72779b8a-f295-4758-b3ec-022b2c20d4e4.png)


## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.


## License
[MIT](https://choosealicense.com/licenses/mit/)
