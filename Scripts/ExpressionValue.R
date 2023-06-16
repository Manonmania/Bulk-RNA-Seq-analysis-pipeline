# Always Load these libraries (Adjust for the final script)
# CRAN
library(assertthat)
library(tidyverse)
library(reshape2)
library(gdata)
library(ggplot2)
library(grid)
library(edgeR)
library(DT)

# Optional

library(knitr)
library(RColorBrewer)
library(gridExtra)
library(flashClust)
library(ggdendro)
library(ggrepel)
library(pheatmap)

project <- readRDS("eset_raw_bulkRNA.RDS")
cell_type <- read.xlsx("Genemarkers.xlsx")

eset_raw <- project$Esets$eset_raw

fullDge <- DGEList(counts = exprs(eset_raw))
fullDge <- calcNormFactors(object = fullDge, method="TMM")

topExprs <- cpm(fullDge, prior.count=3, log=T)

topExprs_Pkp2 <- topExprs[grepl("Pkp2",rownames(topExprs)),] %>% as.data.frame()

colnames(topExprs_Pkp2) <- "Pkp2"

topExprs_Pkp2 <- t(topExprs_Pkp2)
rownames(topExprs_Pkp2) <- "Pkp2"

topExprs_Pkp2_1 <- melt(topExprs_Pkp2)
topExprs_Pkp2_1 <- topExprs_Pkp2_1[,c(2,3)]
colnames(topExprs_Pkp2_1) <- c("Group","Value")

colnames(topExprs_Pkp2_1) <- c("Group","Value")
samples <- colnames(topExprs_Pkp2)
expression <- topExprs_Pkp2[1,]
df <- data.frame(samples,expression)
s <- split(df,df$samples)
Pkp2 <- cbind(s$Vechicle.WT[,2],s$Vehicle.KO[,2])
colnames(Pkp2) <- c("Cre-/- Vehicle","cKO")
write.csv(Pkp2,"Pkp2_expression.csv")