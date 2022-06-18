
source("/home/anna/BIONETS_code/algorithms/GENIE3_R/GENIE3.R") # TODO

args = commandArgs(trailingOnly=TRUE)
data_path = args[0]
regulators = args[1]
targets = args[2]
treeMethod = "RF"
K = "sqrt"
nTrees = args[3]
nCores = 1
returnMatrix = TRUE
verbose = FALSE

exprMatrix <- read.expr.matrix(file_name=data_path, form="rows.are.samples", sep="\t")

weightMat <- GENIE3(
  exprMatrix,
  regulators,
  targets,
  treeMethod,
  K,
  nTrees,
  nCores,
  returnMatrix,
  verbose
)

linkList <- getLinkList(weightMat)
dim(linkList)
head(linkList)

