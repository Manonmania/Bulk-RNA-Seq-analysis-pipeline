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

library(rstatix)
library(ggprism)

df_p_val <- rstatix::t_test(topExprs_Pkp2_1, Value ~ Group, ref.group = "Vehicle.KO") %>% 
  rstatix::add_xy_position()

p <- ggplot(topExprs_Pkp2_1, aes(x = factor(Group), y=Value)) + 
  stat_summary(aes(fill = factor(Group)), geom = "col", fun = mean) +
  stat_summary(geom = "errorbar",
               fun = mean,
               fun.min = function(x) mean(x) - sd(x),
               fun.max = function(x) mean(x) + sd(x),
               width = 0.3) +
  theme_prism() +
  coord_cartesian(ylim = c(0,10)) +
  scale_y_continuous(breaks = seq(0,10,2), expand = c(0,0)) + ggtile("Pkp2") +
  xlab("Group") + ylab("Relative abundance")

# with brackets

p1 <- p + add_pvalue(df_p_val, lable = "adj.p.signif")

p1

# without brackets

p2 <- p + add_pvalue(df_p_val, lable = "p.adj.signif", remove.bracket = TRUE)

p2