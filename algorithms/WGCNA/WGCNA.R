get_netw <- function(prefix) {

        wgcna <- getwd()
        library(WGCNA, warn.conflicts=FALSE)
        enableWGCNAThreads(nThreads=2)
        library(reshape2, warn.conflicts=FALSE)
        main <- file.path(wgcna)

	data_path <- paste(wgcna, '/temp/', prefix, '_expression_data.csv', sep= "")
	out_path <- paste(wgcna, '/temp/', prefix,  '_edge_list.csv', sep = "")

        datExpr <- read.csv(data_path, header = TRUE, sep='\t', as.is=TRUE)
	datExpr <- datExpr[-c(1)]

        # https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/Simulated-05-NetworkConstruction.pdf
        powers = c(c(1:10), seq(from = 10, to=20, by=2))
        # Call the network topology analysis function
        sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
        softPower <- sft$powerEstimate
        if(is.na(softPower)){
            softPower <- 6 # recommended by authors
        }
        line = paste('softPower for ', prefix, ': ', as.character(softPower), ', ', as.character(Sys.time()))
        write(line,file="softPowers.txt",append=TRUE)
        adjacency <-  adjacency(datExpr, power = softPower)

        # Obtaining the edge list from the correlation matrix
        adjacency<-as.matrix(adjacency)
        adjacency[upper.tri(adjacency)]<-NA
        diag(adjacency)<-NA
        edges<-melt(adjacency)
        colnames(edges)<-c("source","target","score")
        edges<-edges[!is.na(edges$score),]

        # Reduce needed storage. Sort out edges with a smaller weight for our purpose, because these won't be in the top 5000 of weighted 
        # edges anyways. This is save, because even if less edges remain, the user will be notified in the saved .csv with the JIs in the
        # column 'filled'. Here, commented.
        #edges<-edges[(edges$score > 0.001),]

        edges$score<-as.numeric(as.character(edges$score))
        edges <- edges[order(edges$score,edges$source,edges$target, decreasing=TRUE),]
	write.table(edges[1:5005,], out_path, sep='\t', row.names=FALSE)
}
                                              
