## Script to estimate association values between nitochondrial and nuclear SNPs using Goodma-Kruskal tau statistic 
### By: Eduardo González
### email: eduardo.glez.or@gmail.com
### INPUT: The input is generated from mitochondrial and nuclear VCF files (1,000 genomes format) using the script genotype_input_generator_GoodMK.NoPhased.pl
### Directory and file names must be change in the script 
library(foreach)
library(data.table)
library(DescTools)
library(parallel)


## Indicate working directory 
setwd("/path/to/working/files/")


### Indicate input files .mtr nuclear (chrmatrix) and mitocondrial (chrMTmatrix)
chrmatrix <- as.matrix(fread("matrix.file.autosomals.mtr", header = FALSE))
chrMTmatrix <- as.matrix(fread("matrix.mitochondrial.mtr", header=FALSE))


### Indicate number of cores to parallel processing 

cl <- makePSOCKcluster(32)



for (x in 1:nrow(chrMTmatrix)) {

	### Asigning names
	Ncol <- ncol(chrmatrix)
	gg1 <- as.vector(chrMTmatrix[x,])
	limit <- Ncol-1
	g1 <- gg1[2:limit]
	nameM <- gg1[1]

	### Exporting variables 

	clusterExport(cl,c('Ncol','g1','GoodmanKruskalTau','nameM'))
	result <- as.matrix(parRapply(cl,chrmatrix,call_GKT))

	#### Output results to .txt

	write.table(result,file=paste(nameM,".outMatrix",sep=","), col.names = FALSE, row.names = FALSE, quote= FALSE)

}



stopCluster(cl)

#####  call_GKT function (.mtr files input needed)

call_GKT <- function(gg2) {
  
  	Gg2 <-as.vector(gg2)
  	limit <- Ncol-1
  	g2 <- Gg2[2:limit]
  	 nameCh <- Gg2[1]
  
  	GG <- GoodmanKruskalTau(g1,g2, direction = c("column"))
    
  	outV <- paste(nameCh,GG,nameM, sep=",")
	   my_mat <- matrix(c(outV), nrow=1,ncol=1)
	  return(my_mat)
}






