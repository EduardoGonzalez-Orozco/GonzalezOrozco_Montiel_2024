# Load required libraries
library(data.table)      # For fast reading and processing of tabular data
library(DescTools)       # For GoodmanKruskalTau and other statistical tools


# Set working directory
setwd("/home/martin/data_lustre/mitoNuclear_diseq/random_test")

# === Load input: allele counts from nuclear and mitochondrial MGF files ===
MGFcountsN <- fread("MXL_all_autosomals_only.snps.mgf.counts.MinVar2.RS")
MGFcountsM <- fread("matrix.MXL.chrMT.phase3_callmom-v0_4.20130502.genotypes.OnlySnps.vcf.NoPhased.mtr.mgf.counts.MinVar2.RS")

# === Initial parameters ===
alleleN <- c("1|1", "0/1")  # Possible genotypes for nuclear SNPs
alleleM <- c("1")           # Possible genotype for mitochondrial SNPs (haploid)

n <- 64        # Number of individuals
s <- 50000000  # Number of simulations
OUT <- rep(".", s)  # Output vector to store Goodman-Kruskal Tau results

# === Create initial vectors with all reference genotypes ===
gN <- rep("0|0", n)  # All nuclear SNPs initially set to homozygous reference
gM <- rep("0", n)    # All mitochondrial SNPs initially set to 0

# === Get allele substitution probabilities from frequency in real data (Nuclear) ===
countS <- as.data.frame(table(MGFcountsN$V6))  # Count occurrences of each substitution level
nonZ <- apply(countS, 1, function(row) all(row != 0))  # Remove rows with 0s
countSClean <- countS[nonZ, ]

# Calculate probabilities
sumFT <- sum(countSClean$Freq)
probs <- countSClean$Freq / sumFT
countSClean[["probs"]] <- probs

SUBSN <- countSClean$Var1   # Substitution numbers
PROBSN <- countSClean$probs # Associated probabilities

# === Repeat for mitochondrial SNPs ===
countSM <- as.data.frame(table(MGFcountsM$V6))
nonZM <- apply(countSM, 1, function(row) all(row != 0))
countSCleanM <- countSM[nonZM, ]

sumFTM <- sum(countSCleanM$Freq)
probsM <- countSCleanM$Freq / sumFTM
countSCleanM[["probs"]] <- probsM

SUBSM <- countSCleanM$Var1
PROBSM <- countSCleanM$probs

# === Main simulation loop ===
startT <- Sys.time()  # Start timing

for (m in 1:s) {
  print(m)

  # Reset genotypes to reference
  gN <- rep("0|0", n)
  gM <- rep("0", n)

  # Sample number of substitutions based on observed probabilities
  SampleNumberN <- sample(SUBSN, 1, prob = PROBSN)
  nx <- round(runif(SampleNumberN, 1, length(gN)))  # Random nuclear positions to change

  SampleNumberM <- sample(SUBSM, 1, prob = PROBSM)
  mx <- round(runif(SampleNumberM, 1, length(gM)))  # Random mitochondrial positions to change

  # Randomly assign new genotypes to sampled positions
  sustN <- sample(alleleN, SampleNumberN, replace = TRUE)
  sustM <- sample(alleleM, SampleNumberM, replace = TRUE)

  # Replace genotypes in vectors
  gN[nx] <- sustN
  gM[mx] <- sustM

  # Compute Goodman-Kruskal's Tau between mitochondrial and nuclear variants
  OUT[m] <- GoodmanKruskalTau(gM, gN, direction = c("column"))
}

# Write output to file
write.table(OUT, file = "MXL.Simulation50M.txt")


