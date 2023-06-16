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

project <- readRDS(file = "eset_raw.RDS")
cell_type <- read.xlsx("Genemarkers.xlsx")

eset_raw <- project$Esets$eset_raw

fullDge <- DGEList(counts = exprs(eset_raw))
fullDge <- calcNormFactors(object = fullDge, method="TMM")

topExprs <- cpm(fullDge, prior.count=3, log=T)

eset_markers <- eset_raw[is.element(fData(eset_raw)$gene_name,
                                    cell_type$Genemarkers), ]

rownames(eset_markers) <- fData(eset_markers)$gene_name

# plotting

anno = merge(fData(eset_markers),
             cell_type,
             by.x = "gene_name",
             by.y = "Genemarkers")

anno <- anno[!duplicated(anno$gene_name),]
rownames(anno) = anno$gene_name
anno <- anno[rownames(eset_markers),]
identical(rownames(anno), rownames(eset_markes))
fData(eset_markers) <- anno

eset_markers <- eset_markers[order(fData(eset_markers)$Populations),]

gaps_samples <- cumsum(as.data.frame(table(paste(pData(eset_markers)$Set)))$Freq)

gaps_genes = cumsum(as.data.frame(table(fData(eset_markers)$Populations))$Freq)

topExprs <- topExprs[rownames(topExprs) %in% rownames(eset_markers),]

reorder_idx <- match(rownames(fData(eset_markers)),rownames(topExprs))
topExprs <- topExprs[reorder_idx,]

gaps_samples <- c(4,10,17)

pheatmap(topExprs,
         scale = "row",
         annotation_row = fData(eset_markers)[, "Populations", drop=FALSE],
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         gaps_col = gaps_samples,
         gaps_row = gaps_genes,
         main = "Heatmap of Genemarkers")