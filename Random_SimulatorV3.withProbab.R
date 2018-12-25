### Código para realizar simulaciones en base a  la variación aleatoria
#### El input se genera a partir del cálculo de la frecuencia mínima alélica por
#### posición en cada población 
#### Escrito por: Eduardo Ma


library(data.table)
library(DescTools)
library(foreach)
library(doParallel)

setwd("/home/martin/data_lustre/mitoNuclear_diseq/random_test")
## MGF file INPUT (counts only)(RS)

MGFcountsN <- fread("MXL_all_autosomals_only.snps.mgf.counts.MinVar2.RS")
MGFcountsM <- fread("matrix.MXL.chrMT.phase3_callmom-v0_4.20130502.genotypes.OnlySnps.vcf.NoPhased.mtr.mgf.counts.MinVar2.RS")
## ParÃ¡metros iniciales

### Alelos posibles
alleleN <-c("1|1","0/1")
alleleM <-c("1")

# Número individuos
n=64

#Número de simulaciones
s=50000000
OUT <- rep(".",s)


## Generate SNPs( full of 0)

gN <- rep("0|0",n)
gM <- rep("0",n)

####### Número de substituciones basado en la frecuencia mínima alélica  (Nuclear):
  
  ##  Obteniendo frecuencia
countS <- as.data.frame(table(MGFcountsN$V6))

  ## Removiendo ceros
nonZ <- apply(countS,1,function(row) all(row !=0))

countSClean <- countS[nonZ,]

  ## Estimating prob
sumFT <- sum(countSClean$Freq)

probs <- (countSClean$Freq/sumFT)

countSClean[["probs"]] <- probs

SUBSN <- countSClean$Var1
PROBSN <- countSClean$probs


####### Número de substituciones basado en la frecuencia mínima alélica  (Mitocondria):

   ##  Obteniendo frecuencia
countSM <- as.data.frame(table(MGFcountsM$V6))

   ## Removiendo ceros

nonZM <- apply(countSM,1,function(row) all(row !=0))

countSCleanM <- countSM[nonZM,]

## Estimating prob
sumFTM <- sum(countSCleanM$Freq)

probsM <- (countSCleanM$Freq/sumFTM)

countSCleanM[["probs"]] <- probsM

SUBSM <- countSCleanM$Var1
PROBSM <- countSCleanM$probs




## Loops  y generación de SNP con variación aleatoria   
startT<- Sys.time()

for (m in 1:s) {
	

   print(m) 
 
  gN <- rep("0|0",n)
  gM <- rep("0",n)
  
  
  # selecting cases
  SampleNumberN <- sample(SUBSN,1,prob = PROBSN)
  nx <- round(runif(SampleNumberN,1,length(gN)))
  
  SampleNumberM <- sample(SUBSM,1,prob = PROBSM)
  mx <- round(runif(SampleNumberM,1,length(gM)))
  
  ## Sampling genotypes to fill (Nuclear)
  sustN <- sample(alleleN,SampleNumberN,replace = TRUE)
  sustM <- sample(alleleM,SampleNumberM,replace = TRUE)
  
  ## Sustituyendo Nuclear
   gN[nx] <-sustN  
  
  ## Loop Mitocondrial
   gM[mx] <-sustM  
  
  ## GoodmanKruskal calculation
  
  OUT[m] <- GoodmanKruskalTau(gM,gN, direction = c("column"))
  
}


write.table(OUT, file ="MXL.Simulation50M.txt")

endT <- Sys.time()

time.taken <- endT- startT

time.taken
