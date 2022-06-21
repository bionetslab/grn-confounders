## Installation instructions

Run these commands in a terminal before running the tests on the GENIE3 tool.
Requires R version 4.2 or older.

```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("GENIE3")
```