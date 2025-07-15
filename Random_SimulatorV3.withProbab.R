# Load required libraries
library(data.table)      # For fast reading and processing of tabular data
library(DescTools)       # For GoodmanKruskalTau and other statistical tools

# Set working directory
setwd("/path/to/working/directory/")

# === Load input: allele counts from nuclear and mitochondrial MGF files ===

# Expected format for MGFcountsN and MGFcountsM:
# These are files generated after processing VCFs into "minor allele count" matrices.
# The files should be tab-delimited, with at least one column containing 
# the number of individuals carrying a minor allele (non-reference) for each SNP.
# Specifically:
# - Column V6 must contain an integer: the number of individuals carrying the minor allele.
# - Each row corresponds to a single SNP.
# - The column V6 (6th column) is the one used here to estimate substitution frequencies.

# Example rows might look like:
# CHR  POS  ...  V6
# 1    1234 ...   4
# 1    2345 ...   1
# 2    3456 ...   5
# (etc.)

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
countS


