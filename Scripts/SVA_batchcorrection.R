# SVA package for removing batch effects and other unwanted variation in high-throughput experiments

# Set the working directory
setwd("/Users/Mano/Documents/CRISPR_Screen_For_AAV5_infection/Downstream_Analysis")

# Load libraries

library(sva)
library(bladderbatch)
library(pamr)
library(limma)
library(dplyr)
library(rafalib)
library(Biobase)
library(DESeq2)
library(edgeR)

# Load the eset
AAA_PROJECT <- readRDS(file = "eset_raw_2023_01_27t00_22_45.RDS")
eset_raw <- AAA_PROJECT$Esets$eset_raw

dge <- DGEList(counts = exprs(eset_raw),
                  remove.zeros = TRUE)

dge <- calcNormFactors(object = dge, method = "TMM")

# Temp Norm for qc
exprs_norm_mat <- cpm(
