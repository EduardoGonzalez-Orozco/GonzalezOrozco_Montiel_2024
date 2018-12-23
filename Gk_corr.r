####Programa para estimar valores de asociaciones entre SNP nucleares y mitocondriales utilizando el estadístico Tau de Goodman-Kruskal 
### Escrito por: Eduardo González
### Correo: Eduardo.glez.or@gmail.com
### Nota: Se necesitan los archivos .mtr de nuclear mitocondria generados por genotype_input_generator_GoodMK.NoPhased.pl

library(foreach)
library(data.table)
library(DescTools)
library(parallel)


## Indicar directorio de trabajo
setwd("/LUSTRE/usuario/martin/mitoNuclear_diseq/CHB")


### Indicar archivos .mtr nuclear (chrmatrix) y mitocondrial (chrMTmatrix)
chrmatrix <- as.matrix(fread(“matrix.file.aotusomals.mtr”, header = FALSE))
chrMTmatrix <- as.matrix(fread("matrix.mitochondrial.mtr", header=FALSE))


### Indicar el número de cores con los que se hará el paralelo el trabajo (En este caso es e 32)

cl <- makePSOCKcluster(32)


### Inicio del ciclo

for(x in 1:nrow(chrMTmatrix) {

	### Asignando nombres a parejas analizar
	Ncol <- ncol(chrmatrix)
	gg1 <- as.vector(chrMTmatrix[x,])
	limit <- Ncol-1
	g1 <- gg1[2:limit]
	nameM <- gg1[1]

	### Exportando variables a trabajar en paralelo 

	clusterExport(cl,c('Ncol','g1','GoodmanKruskalTau','nameM'))
	result <- as.matrix(parRapply(cl,chrmatrix,call_GKT))

	#### Output results to .txt

	write.table(result,file=paste(nameM,".outMatrix",sep=","), col.names = FALSE, row.names = FALSE, quote= FALSE)

}



stopCluster(cl)

##### Función call_GKT la cual. Recibe como input la matriz (archivos .mtr)

call_GKT <- function(gg2) {
  
  	Gg2 <-as.vector(gg2)
  	limit <- Ncol-1
  	g2 <- Gg2[2:limit]
  	 nameCh <- Gg2[1]
  
  	GG <- GoodmanKruskalTau(g1,g2, direction = c("column"))
    
  	outV <- paste(nameCh,GG,pvalue,nameM, sep=",")
	   my_mat <- matrix(c(outV), nrow=1,ncol=1)
	  return(my_mat)
}






