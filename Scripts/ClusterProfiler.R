# Set current directory

setwd("/home/mano/project")

library(clusterProfiler)
library(pathview)
library(enrichplot)
library(ggplot2)
library(DESeq2)

organism = "org.Mm.eg.db"
library(organism, character.only = TRUE)

eset <- readRDS("./save/Preliminary_contrast.RDS")

topTableFull <- eset$Analysis[["Group Analysis"]]$DEG$results_full
KO_WT <- topTableFull[,13:14]
KO_WT_Dn <- KO_WT[KO_WT$`Vehicle.KO vs. Vehicle.WT_Regulation` == "Dn",]
sym <- rownames(KO_WT_Dn)
EG_IDs = mget(sym, revmap(org.Mm.egSYMBOL),ifnotfound = NA)
KEGG_IDs = mget(as.character(EG_IDs), org.Mm.egPATH,ifnotfound = NA)

original_gene_list <- KO_WT_Dn$`Vehicle.KO vs. Vehicle.WT_logFC`
names(original_gene_list) <- rownames(KO_WT_Dn)

ids <- bitr(names(original_gene_list), fromType = "SYMBOL",
            toType = "ENTREZID", OrgDb=organism)
dedup_ids <- ids[!duplicated(ids[c("SYMBOL")]),]
df2 = KO_WT_Dn[rownames(KO_WT_Dn) %in% dedup_ids$SYMBOL,]
df2$Y <- dedup_ids$ENTREZID
kegg_gene_list <- df2$`Vehicle.KO vs. Vehicle.WT_logFC`
names(kegg_gene_list) <- df2$Y
kegg_gene_list <- na.omit(kegg_gene_list)
kegg_gene_list = sort(kegg_gene_list, decreasing =TRUE)
kegg_organism = "mmu"

kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = kegg_organism,
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")

kk2@result$Description <- gsub("Mus musculus","",kk2@result$Description)
kk2@result$Description <- gsub("\\(house mouse)","",kk2@result$Description)
kk2@result$Description <- gsub("-","",kk2@result$Description)

library(ggplot2)
library(dplyr)
library(stringr)

# count the gene number

gene_count <- kk2@result %>% group_by(ID) %>% summarise(count = sum(str_count(core_enrichment, "/")) +1)

# merge with the original dataframe

dot_df <- left_join(kk2@result, gene_count, by = "ID") %>% mutate(GeneRatio = count/setSize)
dot_df1 <- dot_df[dot_df$NES < 0,]

# plot
library(forcats)
ggplot(dot_df1, aes(x = GeneRatio, y = fct_reorder(Description, GeneRation))) +
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  theme_bw(base_size = 14) +
  scale_colour_gradient(limits=c(0,0.10), low="res") +
  ylab(NULL) +
  ggtitle("Vehicle WT vs. Vehicle KO")

dotplot(kk2, showCategory = 10, title = "Vehicle WT vs. Vehicle KO", split = ".sign") +
  facet_grid(.~.sign)
