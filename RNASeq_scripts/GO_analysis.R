################# GO Analysis

## Load data from DEGs analysis
# Read data
DEGs = "RUV_Normalization_27042016/HSC_DEGs_annotation.xls"
all_gene = "RUV_Normalization_27042016/HSC_DESeq_Results.xls"

# DEGs
df = read.csv(DEGs, sep = "\t", header = T)
# Extract IDs from differential expression analysis
ensID = unique(df[,1])
ensID_UP = unique(df[which(df$log2FoldChange>0),1])
ensID_DOWN = unique(df[which(df$log2FoldChange<0),1])

eg = unique(df$entrezgene[!is.na(df$entrezgene)])
eg_UP = unique(df[!is.na(which(df$log2FoldChange>0)),"entrezgene"])
eg_DOWN = unique(df[!is.na(which(df$log2FoldChange<0)),"entrezgene"])

# All genes
df_all = read.csv(all_gene, sep = "\t", header = T)
# Translate ENS_ID in ENTREZ
all_ens = df_all[,1]
all_eg = bitr(all_ens, fromType="ENSEMBL", toType="ENTREZID", annoDb="org.Mm.eg.db")
all_log2FC = merge(all_eg, df_all, by.x="ENSEMBL", by.y="ensembl_gene_id")[,c("ENTREZID","log2FoldChange")]
val=as.vector(all_log2FC$log2FoldChange)
names(val)=all_log2FC$ENTREZID
## val = geneList for KEGG and GSEA analysis
val=sort(val, decreasing = T)
head(val)
### Different methods

# 1 - clusterProfiler
# Following bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html

# Load libraries
library(clusterProfiler)

## Gene Ontology Classification
library("DOSE")
# Select the gene list to use for enrichment
gene = eg_DOWN
geneList = val
ont = "BP"

# Genes grouping by GO
ggo <- groupGO(gene     = as.character(gene),
               organism = "mouse",
               ont      = ont,
               level    = 4,
               readable = TRUE)
head(summary(ggo))
barplot(ggo)
# GO over-representation test
ego <- enrichGO(gene          = gene,
                organism      = "mouse",
                ont           = ont,
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.001,
                readable      = TRUE)
#head(summary(ego))
#dotplot(ego)
res=ego@result
nrow(res)
# KEGG Analysis
kk <- enrichKEGG(gene         = gene[which(!is.na(as.character(gene)))],
                 organism     = "mouse",
                 pvalueCutoff = 0.05, 
                 readable     = TRUE,
                 use_internal_data = TRUE)
head(summary(kk))

# DAVID enrichment analysis
david <- enrichDAVID(gene = gene,
                     idType = "ENTREZ_GENE_ID",
                     listType = "Gene",
                     annotation = "KEGG_PATHWAY",
                     david.user = "clusterProfiler@hku.hk")

# KEGG Gene set analysis
kk2 <- gseKEGG(geneList     = geneList,
               organism     = "mouse",
               nPerm        = 1000,
               minGSSize    = 120,
               pvalueCutoff = 0.01,
               verbose      = FALSE,
               use_internal_data = TRUE)


# GO Semantic similarity
library(GOSemSim)
library(pheatmap)
go1= rownames(res)
go2= rownames(res)

goMatrix=mgoSim(go1, go2, ont=ont, combine=NULL)
# Removing ROOT MF
goMatrix=goMatrix[!rownames(goMatrix)=="GO:0003674",!rownames(goMatrix)=="GO:0003674"]
#create the breaks
bk2 = unique(c(seq(0, 0.39, length=5), 0.5, seq(0.51,1, length=15)))
#set different color vectors for each interval
col1 = colorRampPalette(c("white", 'orange'))(5) #set the order of greys
col2 <- rep("orange", 1)
col3 = colorRampPalette(c("orange", "red"))(15)
colors2 <- c(col1, col2, col3)

pheatmap(goMatrix, breaks = bk2 , color = colors2)

a=pheatmap(goMatrix, breaks = bk2 , color = colors2, kmeans_k = 2)
write.table(names(which(a$kmeans[[1]]==1)),"test", quote = F, row.names = F, col.names = F)

hc=hclust(dist(goMatrix))