## GRN confounder -- draft --
Assessing effect of sex, age, and ethnicity confounders on the GRN inferred by the [GENIE3 tool](https://github.com/vahuynh/GENIE3).

GENIE3 is a tool for GRN inference that uses tree-based ensemble methods. The paper in which the tool was introduced first by the developers ... can be found [here]( 10.1371/journal.pone.0012776).
The subproject located in this directory of the repository aims to assess the robustness of GENIE3 with respect to data biases induced by the confounders sex, age and ethnicity.

Notes about what is important to mention here later:
- parameters were selected based on the recommendations in the paper/according to the parameters that produced the best results (benchmark: DREAM4 in-silico data set)
- difference between race and ethnicity in TCGA and the reason for that particular separation
- arguments one has to hand to script
- input/output (format; what does preprocessing do)
- most important steps in script: only Primary Tumor, get confounder based partition on samples, get random partitions
- run GENIE3 multiple times also on the permanent confounder based partition, becuse GENIE3 includes randomized procedures
- run GENIE3 on n random partitions
- normalize and remove 0-cols and explain why; sync genes in data, gene list and regulators

usage:
Call the script data_partitioning.py. When calling the script for the first time on a new system, use the option -prp to enable data preprocessing. In all subsequent executions of the script, the -prp option can be omitted as the preprocessed data is available locally in .pkl format.
![usage](https://user-images.githubusercontent.com/62886428/172936096-da788d0d-3d2c-4412-ba18-765594d6b4e1.png)

Then explain briefly data visualization
