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
library(EnhacedVolcano)

project <- readRDS("eset_raw_bulkRNA.RDS")
eset_raw <- project$Esets$eset_raw

lst_contrast <- readRDS("save/Preliminary_contrast.RDS")

theContrast <- lst_contrast$Analysi$`Group Analysis`
fullEset <- theContrast$DEG$eset
topExprs <- exprs(eset_raw)
topTableFull <- theContrast$DEG$result_full
gene.list <- c()

topTableTemp <- topTableFull[,c(12,14)] %>% as.data.frame()
colnames(topTableTemp) <- c("adj.P.Val","log2FoldChange")

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'adj.P.Val',
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 10e-32,
                FCcutoff = 2.0,
                pointSize = 4.0,
                labSize = 6.0,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75)

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'adj.P.Val',
                selectLab = c("Jup","Pkp2")
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 10e-32,
                FCcutoff = 2.0,
                pointSize = 4.0,
                labSize = 6.0,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75)