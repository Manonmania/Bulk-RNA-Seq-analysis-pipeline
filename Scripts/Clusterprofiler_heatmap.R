
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

result <- kk2@result %>% as.data.frame()
rownames(result) <- NULL

gene_name <- gsub("/"," ",result$core_enrichment)
df2$gene_name <- rownames(df2)

gene_set <- list()
for (j in gene_name){
  indi_gene <- strsplit(j," ") %>% as.data.frame()
  for (k in indi_gene){
    gene_res <- df2[df2$Y %in% k,]
  }
  gene_set[[j]] <- gene_set$gene_name
}

result_gene <- list()
for (l in 1:length(gene_set)){
  result_gene[[l]] <- toString(gene_set[[l]])
}

result_gene <- unlist(result_gene)
library(openxlsx)
cell_type <- read.xlsx("ARVC_genemarkers.xlsx")

project <- readRDS("eset_raw_bulkRNA.RDS")
eset_raw <- project$Esets$eset_raw

library(edgeR)

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

topExprs <- topExprs[rownames(topExprs) %in% rownames(eset_markers),]

library(pheatmap)
pheatmap(topExrs, cluster_rows = FALSE, cluster_cols = FALSE, scale = "row")