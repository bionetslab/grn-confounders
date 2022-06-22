## Installation instructions

Run these commands in a terminal before running the tests on the GENIE3 tool.
Requires R version 4.2 or older.

```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("GENIE3")
```

or try:
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("GENIE3")
```

[Reference](https://bioconductor.org/packages/release/bioc/html/GENIE3.html)

Note: the message "The following object is masked from  .GlobalEnv: GENIE3" refers to a testing function in GENIE3 that is not required here. 

