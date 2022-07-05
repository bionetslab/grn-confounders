
init = function(workingDir) {

getwd();
setwd(workingDir);

library(WGCNA);
library(tidyverse)
library(xlsx)
library(openxlsx);
library("readxl")
library(rlist)
library(dplyr)
options(stringsAsFactors = FALSE);

enableWGCNAThreads()
}

preProcess = function(geneFile, phenotypeFile, confounder, nRandomPartitions){

geneData <<- read.csv(file = geneFile, head = TRUE, sep="\t");
phenotypeData <<- read.csv(file = phenotypeFile, head = TRUE, sep="\t")
confounder <<- confounder;
nRandomPartitions <<- nRandomPartitions;

if(confounder == "race"){
	confounderData <<- split(phenotypeData, phenotypeData$race.demographic)
}

if(confounder == "ethnicity"){
	confounderData <<- split(phenotypeData, phenotypeData$ethnicity.demographic)
}

if(confounder == "sex"){
	confounderData <<- split(phenotypeData, phenotypeData$gender.demographic)
}
}

###

GeneDataCleanUp = function(geneData){

n <- geneData[,1]
transposedGeneData <<- as.data.frame(as.matrix(t(geneData[,-1])))
colnames(transposedGeneData)<- n

transposedGeneData <<- cbind(rownames(transposedGeneData), transposedGeneData)
rownames(transposedGeneData) <<- NULL
colnames(transposedGeneData)[1] <<- "submitter_id.samples"
#transposedPancreasCancerData[c(1:5),c(1:5)]
#dim(transposedPancreasCancerData)
}
###

###

confounderGeneDataCleanUp = function(confounderData, transposedGeneData){

Jaccard_Indexes_cnf <<- data.frame(matrix(ncol = 2, nrow = 0))
x <<- c("Confounder_partition", "Jaccard_Index")
colnames(Jaccard_Indexes_cnf) <<- x
ConfounderGenesList <<- data.frame(matrix(ncol = 1, nrow = 0))
y <<- c("submitter_id.samples")
colnames(ConfounderGenesList) <<- y
confoundersAdjacencyList <<- vector(mode = "list", length = 0)
confoundersLength <<- vector(mode = "list", length = 0)
for (i in 1:length(confounderData)) {
  if(names(confounderData[i]) != "not reported"){
	confounderData[[i]] = confounderData[[i]][c("submitter_id.samples")]
	#print(dim(confounderData[[i]]))
	#print(length(confoundersLength))
	confounderData[[i]] = data.frame(lapply(confounderData[[i]], function(x) {gsub("-", ".", x)}))
	#print(confounderData[[i]])
	mergedData = merge(confounderData[[i]], transposedGeneData, by =c("submitter_id.samples"))
	ConfounderGenesList <<- rbind(ConfounderGenesList ,mergedData[c(1)])
	confoundersLength <<- list.append(confoundersLength, dim(mergedData)[1])

	mergedData = mergedData[-c(1)]
	#print(dim(mergedData))

	confAdjacency <<- findAdjacency(mergedData[c(1:dim(mergedData)[1]),c(1:20)])
	edgesFirstBlock <<- edgesFromAdjacency(confAdjacency)

	confoundersAdjacencyList <<- list.append(confoundersAdjacencyList, edgesFirstBlock)
}
}

	if(length(confoundersAdjacencyList) == 2){
	findJaccard = findJaccardIndex2Matrix(confoundersAdjacencyList[[1]], confoundersAdjacencyList[[2]])
	print(findJaccard)
	Jaccard_Indexes_cnf[nrow(Jaccard_Indexes_cnf) + 1,] <<- c(confounder, findJaccard)
	}

	if(length(confoundersAdjacencyList) == 3){
	findJaccard = findJaccardIndex2Matrix(confoundersAdjacencyList[1], confoundersAdjacencyList[2], confoundersAdjacencyList[3])
	print(findJaccard)
	#Jaccard_Indexes_cnf[nrow(Jaccard_Indexes_cnf) + 1,] <<- c(cnf, findJaccard)
	}
}

###


findAdjacency = function(Genedata){

Corr <- adjacency(Genedata, type = "signed", power = 1)
#Corr[c(1:8),c(1:8)]
Adjacency <- signumAdjacencyFunction(Corr, threshold = 0.6)

return(Adjacency)

}

###

randomPartitions = function(nRandomPartitions, confoundersLength, phenotypeData, transposedGeneData){

Jaccard_Indexes_rnd <<- data.frame(matrix(ncol = 2, nrow = 0))
x <<- c("Rnd_Partition_Num", "Jaccard_Index")
colnames(Jaccard_Indexes_rnd) <<- x

phenotypeData <<- phenotypeData
#phenotype <<- phenotypeData[c("submitter_id.samples")]
#phenotype <<- data.frame(lapply(phenotype, function(x) {gsub("-", ".", x)}))
#print(dim(phenotype))
if (length(confoundersLength) == 2){
	n1 <<- confoundersLength[[1]]
	n2 <<- confoundersLength[[2]]
for (i in 1:nRandomPartitions) {
	#RandomisedPhenotype <<- setNames(data.frame(matrix(ncol = 1, nrow = (n1+n2))), c("submitter_id.samples"))
	#rand = phenotype[sample(nrow(phenotypeData)),]
	RandomisedPhenotype <<- as.data.frame(matrix(phenotypeData[sample(1:nrow(phenotypeData)),]), ncol= 1, byrow = TRUE)
	names(RandomisedPhenotype)[1] <<- "submitter_id.samples"
	firstBlock <<- setNames(data.frame(matrix(ncol = 1, nrow = n1)), c("submitter_id.samples"))
	secondBlock <<- setNames(data.frame(matrix(ncol = 1, nrow = n2)), c("submitter_id.samples"))
	firstBlock$submitter_id.samples <<-  RandomisedPhenotype[1:n1,]
	secondBlock$submitter_id.samples <<- RandomisedPhenotype[(n1+1):(n1+n2),]
	firstBlockGenes <<- merge(firstBlock , transposedGeneData, by =c("submitter_id.samples"))
	firstBlockGenes <<- firstBlockGenes[-c(1)]
	secondBlockGenes <<- merge(secondBlock, transposedGeneData, by =c("submitter_id.samples"))
	secondBlockGenes <<- secondBlockGenes[-c(1)]

	firstBlockAdjacency <<- findAdjacency(firstBlockGenes[c(1:dim(firstBlockGenes)[1]),c(1:20)])
	edgesFirstBlock <<- edgesFromAdjacency(firstBlockAdjacency)
	secondBlockAdjacency <<- findAdjacency(secondBlockGenes[c(1:dim(secondBlockGenes)[1]),c(1:20)])
	edgesSecondBlock <<- edgesFromAdjacency(secondBlockAdjacency)

	print(edgesFirstBlock)
	print(edgesSecondBlock)

  	findJaccard = findJaccardIndex2Matrix(edgesFirstBlock, edgesSecondBlock)
	Jaccard_Indexes_rnd[nrow(Jaccard_Indexes_rnd) + 1,] <<- c(i, findJaccard)
}
}
else {
	n1 = confoundersLength[[1]]
	n2 = confoundersLength[[2]]
	n3 = confoundersLength[[3]]

for (i in 1:nRandomPartitions) {

	RandomisedPhenotype <<- as.data.frame(matrix(phenotypeData[sample(1:nrow(phenotypeData)),]), ncol= 1, byrow = TRUE)
	names(RandomisedPhenotype)[1] <<- "submitter_id.samples"

	firstBlock <<- setNames(data.frame(matrix(ncol = 1, nrow = n1)), c("submitter_id.samples"))
	secondBlock <<- setNames(data.frame(matrix(ncol = 1, nrow = n2)), c("submitter_id.samples"))
	thirdBlock <<- setNames(data.frame(matrix(ncol = 1, nrow = n2)), c("submitter_id.samples"))

	firstBlock$submitter_id.samples <<-  RandomisedPhenotype[1:n1,]
	secondBlock$submitter_id.samples <<- RandomisedPhenotype[(n1+1):(n1+n2),]
	thirdBlock$submitter_id.samples <<- RandomisedPhenotype[(n1+n2+1):(n1+n2+n3),]
	firstBlockGenes <<- merge(firstBlock , transposedGeneData, by =c("submitter_id.samples"))
	firstBlockGenes <<- firstBlockGenes[-c(1)]
	secondBlockGenes <<- merge(secondBlock, transposedGeneData, by =c("submitter_id.samples"))
	secondBlockGenes <<- secondBlockGenes[-c(1)]
	thirdBlockGenes <<- merge(thirdBlock, transposedGeneData, by =c("submitter_id.samples"))
	thirdBlockGenes <<- thirdBlockGenes[-c(1)]


	firstBlockAdjacency <<- findAdjacency(firstBlockGenes[c(1:dim(firstBlockGenes)[1]),c(1:20)])
	edgesFirstBlock <<- edgesFromAdjacency(firstBlockAdjacency)

	secondBlockAdjacency <<- findAdjacency(secondBlockGenes[c(1:dim(secondBlockGenes)[1]),c(1:20)])
	edgesSecondBlock <<- edgesFromAdjacency(secondBlockAdjacency)

	thirdBlockAdjacency <<- findAdjacency(thirdBlockGenes[c(1:dim(thirdBlockGenes)[1]),c(1:20)])
	edgesThirdBlock <<- edgesFromAdjacency(thirdBlockAdjacency)

	findJaccard = findJaccardIndex3Matrix(edgesFirstBlock, edgesSecondBlock, edgesThirdBlock)
	Jaccard_Indexes_rnd[nrow(Jaccard_Indexes_rnd) + 1,] <<- c(i, findJaccard)


}
}
}
###


edgesFromAdjacency = function(AdjacencyMatrix){
Edges<<-graph_from_adjacency_matrix(AdjacencyMatrix, mode = c("undirected"))
Edges<<- as.data.frame(get.edgelist(Edges))
Edges$edges <<- paste(Edges[,c(1)], Edges[,c(2)], sep = "-")
Edges <<- Edges[c(3)]

return(Edges)
}

findJaccardIndex2Matrix = function(Edges_A, Edges_B){

intersection = dim(intersect(Edges_A, Edges_B))[1]
union = dim(Edges_A)[1] + dim(Edges_B)[1] - intersection
return (intersection/union)

#Graph_A<<-graph_from_adjacency_matrix(Adjacency_A, mode = c("undirected"))
#Graph_B<<-graph_from_adjacency_matrix(Adjacency_B, mode = c("undirected"))
#Graph_A <<- get.adjacency(Graph_A)
#Graph_B <<- get.adjacency(Graph_B)
#A<<-sum(Graph_A != Graph_B)
#B<<-sum(Graph_A * Graph_B)
#return(round(B/sum(A,B),digits = 4))

}

findJaccardIndex3Matrix = function(Edges_A, Edges_B, Edges_C){


intersection = dim(intersect(Edges_A, Edges_B, Edges_C))[1]
union = dim(Edges_A)[1] + dim(Edges_B)[1] + dim(Edges_C)[1] - intersection
return (intersection/union)

}

#######



###-MAIN-###

init("C:/Users/Mohamed Ughratdar/OneDrive/Desktop/FAU/Project-DataBias/PancreaticCancerDataAnalysis")

preProcess("TCGA-PAAD.htseq_fpkm.tsv", "TCGA-PAAD.GDC_phenotype.tsv", "sex", 2)

GeneDataCleanUp(geneData)

confounderGeneDataCleanUp(confounderData, transposedGeneData)

randomPartitions(1, confoundersLength, ConfounderGenesList, transposedGeneData)
###

