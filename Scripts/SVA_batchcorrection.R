# SVA package for removing batch effects and other unwanted variation in high-througput experiments

# Set current directory

setwd("/home/mano/project")

# load libraries

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

project <- readRDS("eset_raw_bulkRNA.RDS")
eset_raw <- project$Esets$eset_raw

dge <- DGEList(counts = exprs(eset_raw),
               remove.zeros = TRUE)
dge <- calcNormFactors(object = dge, method="TMM")

# Temp Norm for qc

exprs_norm_mat <- cpm(dge[rowSums(dge$counts) > dim(dge)[2],], normalized.lib.sizes = TRUE, log = T, prior.count = 0.25)
eset_filter <- eset_raw[is.element(fData(eset_raw)$gene_name, rownames(exprs_norm_mat)),]

eset_norm <- ExpressionSet(assayData = as.matrix(exprs_norm_mat),
                           phenoData = new("AnnotatedDataFrame", data=pData(eset_raw)),
                           featureData = new("AnnotatedDataFrame", data=fData(eset_filter)),
                           annotation = "Homo")

# Construct DESeqDataSet Object

dds <- DESeqDataSetFrameMatrix(countData = exprs(eset_raw),
                               design ~ 1,
                               colData = pData(eset_raw))
dds <- dds[rowMeans(counts(dds, normalized=FALSE)) > 5,]

# vst function perform variance stabilizing transformation

vsd = vst(dds, blind = TRUE, nsub = 10000)
norm_counts <- assay(vsd)

# Remove unwanted columns in the sample data information
pdata <- pData(eset_raw)
pdata1 <- pdata[,c(2,4,5,11)]

# Only uninfected samples

pdata1$treatment <- ifelse(pdata1$Group == "S","Control","KO_genes")

pdata1$sample_name <- gsub("_Rep[0-9","",pdata1$sample_name)

# Create the full model matrix of variable of interest
mod = model.matrix(~sample_name, data=pdata1)

colnames(mod) <- gsub(pattern = "sample_name",
                      replacement = "",
                      colnames(mod))

# The null model contains only the adjustment variables. Sincer we are not adjusting for any other
# variables in this analysis, only an intercept is included in the model

mod0 = model.matrix(~1, data=pdata1)

# Apply the sva function to estimate the surrgoate variables
svobj=sva(exprs_norm_mat,mod,mod0)
#svobj = sva(norm_counts,mod,mod0)
full.model.sv <- cbind(mod, svobj$sv)

# Fit the linear model with the surrogate variables included
fit <- lmFit(exprs_norm_mat, full.model.sv)

# Get clean data for PCA:

# Regress out all surrgoate variables and keep variance caused by group 
# which are columns 2 to 48, and intercept, which is column 1):

mod_PCA <- coefficients(fit)[,-c(1:48)] %*% t(fit$design[,-c(1:48)])
eset_PCA <- eset_norm
exprs(eset_PCA) <- exprs_norm_mat - mod_PCA

pca_res <- prcomp(t(exprs(eset_PCA)), scale = TRUE)
s <- summary(pca_res)$important[, 1:5]
pca_scores <- pca_res$x %>% as.data.frame()
var_explained <- pca_res$sdev^2/sum(pca_res$sdev^2)

# PCA plot

library(ggplot2)
library(ggrepel)

plotDF <- data.frame(Dim1 = pca_scores$PC1, Dim2 = pca_score$PC2,
                     Group = eset_raw$treatment,
                     Lable = eset_raw$name)

ggplot(data = plotDF,
       mapping = aes(x = Dim1,
                     y = Dim2,
                     color = Group,
                     lable = Label)) +
  theme_bw(base_size=32) +
  labs(x=paste0("PC1: ", round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ", round(var_explained[1]*100,1),"%")) +
  geom_point(size = 5) + geom_repel()

saveRDS(eset_PCA, file=paste0("eset_batcheffects_removed",".RDS"))
