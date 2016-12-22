#################
####         ####
#### DESEQ2  ####
####         ####
#################

library("BiocParallel")
register(MulticoreParam(10))
library('DESeq2')

csv_path='SampleSheet_HSC_outlierOut.tsv'
ref_condition="WT"

### Example sample sheet
# Sample sheet is SPACE-DELIMITED, with header
# Example
# sampleName fileName condition1
# HET_721 HET_721.HQ.counts HET
# HET_729 HET_729.HQ.counts HET

sampleTable <- read.table(csv_path,header=TRUE)
sample_dir='../mapping/'

## DESIGN
# Single factor -> CONDITION
#ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,directory = sample_dir, design= ~ condition)
# MULTI FACTOR -> Condition + proliferation status
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,directory = sample_dir, design= ~ LogCellCount + condition)
##

# DESeq worker
dds <- DESeq(ddsHTSeq, parallel=TRUE)
# Plot dispersion
plotDispEsts(dds)

## Pre-filtering
# All samples > 1 in gene count
dds<-dds[rowSums(counts(dds)) > 1,] 
# Extracting counts
raw.counts <- as.data.frame(counts( dds ))
normalized.counts <- as.data.frame(counts( dds, normalized=TRUE ))

# Relevel the condition
dds$condition <- relevel(dds$condition, ref=ref_condition)
res <- results(dds, parallel=TRUE)
head(res)
# Summary of DEGs
summary(res)


#stabilizing the variance (takes long...)
rld <- rlog(dds)


# Variable for the condition
populations = levels(dds$condition)

#Generate 'means' dataframe
len = length(rowMeans(counts(dds,normalized=TRUE)[,dds$condition == populations[19]]))
means = data.frame(matrix(NA, nrow = len))
count = 0
for (pop in populations){
  print(pop)
  tryCatch({means[pop] = rowMeans(counts(dds,normalized=TRUE)[,dds$condition == pop])},
           error = function(e) {
             means[pop] = counts(dds,normalized=TRUE)[,dds$condition == pop]
           })    
}
# Remove NA initial column from 'means'
means = means[,2:length(means)]
rownames(means) = rownames(normalized.counts)

means$ensembl_gene_id = rownames(means)

#retrieve gene data
library("biomaRt")
gene_set = rownames(means)
ensembl=useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "mmusculus_gene_ensembl",
                host = "feb2014.archive.ensembl.org",
                path = "/biomart/martservice")
genes_information = getBM(c("mgi_symbol",
                            "entrezgene",
                            "uniprot_genename",
                            "description",
                            "gene_biotype",
                            "chromosome_name",
                            "strand",
                            "start_position",
                            "end_position",
                            "ensembl_gene_id"),
                          "ensembl_gene_id",
                          gene_set, ensembl)
genes_with_information = merge(x = means, y = genes_information, by = "ensembl_gene_id", all.x = TRUE)

# Select significative genes
alpha=0.1
resSig <- res[which(res$padj <=alpha),]
downDEGs <- res[which(res$padj <=alpha & res$log2FoldChange < 0),]
upDEGs <- res[which(res$padj <=alpha & res$log2FoldChange > 0),]
# Which genes are significant
ens_genesSig=rownames(resSig)
# Extract annotation
geneSig = genes_with_information[genes_with_information$ensembl_gene_id %in% ens_genesSig,]
# Create DF for means with geneSig information
means$significant[means$ensembl_gene_id %in% ens_genesSig] <- 1
means$significant[!(means$ensembl_gene_id %in% ens_genesSig)] <- 0
genes_with_information$significant[genes_with_information$ensembl_gene_id %in% ens_genesSig] <- 1
genes_with_information$significant[!(genes_with_information$ensembl_gene_id %in% ens_genesSig)] <- 0

# Write DEGs with log2FC
#DEG_list = cbind(as.data.frame(downDEGs[2]),ensid=rownames(downDEGs))
#merge(DEG_list,geneSig,by.x="ensid",by.y="ensembl_gene_id")
#write.table(rownames(downDEGs),"test.rnk", quote=F ,row.names = F, col.names = F, sep="\t")

###VISUALIZATION

# Plot contus of specific genes
library(ggplot2)
### List of specific genes
## From Marica presentation
gene_list = c("TNF","IL1B","IL6","CXCL1","CXCL10","INSR","LOXL2","CCL2","COL1A1","ACTA2","INFG")

ensgene_list = sapply(gene_list, function(x){
  genes_with_information[which(genes_with_information$uniprot_genename==x),]$ensembl_gene_id
  })
# Plot gene counts
par(mfrow=c(3,1),mar=c(3,5,3,1))
sapply(names(ensgene_list),function(x){
  if(length(ensgene_list[[x]]) != 0){
    plotCounts(dds,ensgene_list[[x]], main=paste(as.character(x),as.character(ensgene_list[[x]]), sep = " - "))
  }
})

# ggplot version
#d = plotCounts(dds,"ENSMUSG00000005534", returnData = TRUE)
#ggplot(d, aes(x=condition,y=count)) +
#  geom_point(position = position_jitter(w=0.1,h=0))+
#  ggtitle("INSR")+
#  theme(plot.title = element_text(lineheight=.8, face="bold"))


## Following DESeq 2 tutorial
sampleDists <- dist( t( assay(rld) ) )
## Clssical sample distance
library("pheatmap")
library("gplots")
library("RColorBrewer")
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <-  colnames(rld)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

#sample distances with names
sampleDists <- dist( t( assay(rld) ) )
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(colnames(sampleDistMatrix), sep="-" )
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
hc <- hclust(sampleDists)
heatmap.2(sampleDistMatrix, Rowv=as.dendrogram(hc),
           symm=TRUE, trace="none", col=colors,
           margins=c(2,10), labCol=TRUE )

# Sample distances with Poisson distance calculation
library("PoiClaClu")
poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
rownames(samplePoisDistMatrix) <- paste( colnames(rld), rld$LogCellCount, sep="-" )
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows=poisd$dd,
         clustering_distance_cols=poisd$dd,
         col=colors)

#sample PCA
plotPCA(rld, intgroup = "condition")

# Classic distance
data_PCA <- plotPCA(rld, intgroup = c("condition","LogCellCount"), returnData=TRUE)
percentVar <- round(100 * attr(data_PCA, "percentVar"))
library("ggplot2")
qplot(PC1, PC2, color=LogCellCount, data=data_PCA, label=rownames(data_PCA)) +
  scale_colour_gradient(low="black", high="red") +
  geom_point(aes (shape = factor(condition) ), size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + geom_text(vjust = 1.7, nudge_x = 0.05, check_overlap = T)

## Scatter plot for mean values
library(ggplot2)
colScale <- scale_colour_manual(name = "significant",values = c("gray","red"))
ggplot(genes_with_information, aes(WT, HET, colour=as.factor(significant), alpha = as.factor(significant))) +
  geom_point() +
  colScale +
  # label
  #geom_text(aes(label=ifelse((significant==1),toupper(as.character(mgi_symbol)),'')),
  #          hjust=0, 
  #          vjust=0, 
  #          check_overlap = T) +
  scale_y_log10() +
  scale_x_log10()

## Gene clustering
# Top Variable genes
library("genefilter")
topVarGenes <- head(order(rowVars(assay(rldssva)),decreasing=TRUE),30)
mat <- assay(rldssva)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
# Create Gene name dataframe
df_gene <- as.data.frame(rownames(mat))
colnames(df_gene) <- "ensembl_gene_id"
df_gene = merge(x = df_gene, y = genes_information, by = "ensembl_gene_id", all.x = TRUE)
# fix NA of BLANK values with ENSID
idx_to_fix = which(df_gene$uniprot_genename == "" | is.na(df_gene$uniprot_genename))
df_gene[idx_to_fix,]$uniprot_genename = as.vector(df_gene[idx_to_fix,]$ensembl_gene_id)
rownames(mat) = df_gene[match(df_gene$ensembl_gene_id,rownames(mat)),]$uniprot_genename
# Create condition dataframe
df <- as.data.frame(colData(rldssva)[,c("type","condition")])
# Plot the heatmap
pheatmap(mat, annotation_col=df)

# DEGs Gene
alpha=0.1
resSig <- res[which(res$padj <=alpha),]

mat <- assay(rld)[ which(rownames(rld) %in% rownames(resSig)),]
mat <- mat - rowMeans(mat)
# Create Gene name dataframe
df_gene <- as.data.frame(rownames(mat))
colnames(df_gene) <- "ensembl_gene_id"
df_gene = merge(x = df_gene, y = genes_information, by = "ensembl_gene_id", all.x = TRUE)
# fix NA of BLANK values with ENSID
idx_to_fix = which(df_gene$mgi_symbol == "" | is.na(df_gene$mgi_symbol))
df_gene[idx_to_fix,]$mgi_symbol = as.vector(df_gene[idx_to_fix,]$ensembl_gene_id)
rownames(mat) = unique(df_gene[match(df_gene$ensembl_gene_id,rownames(mat)),]$mgi_symbol)
# Create condition dataframe
df <- as.data.frame(colData(rld)[,"condition"])
colnames(df) <- "condition"
rownames(df) <- rownames(colData(rld))
# Plot the heatmap
pheatmap(mat, annotation_col = df)


### REMOVING Batch effect
## Following DESeq2 tutorial
library("sva")
dat <- counts(dds, normalized=TRUE)
idx <- rowMeans(dat) > 1
dat <- dat[idx,]
mod <- model.matrix(~ condition, colData(dds))
mod0 <- model.matrix(~ 1, colData(dds))
svseq <- svaseq(dat, mod, mod0, n.sv=2)

svseq$sv

par(mfrow=c(2,1),mar=c(3,5,3,1))
stripchart(svseq$sv[,1] ~ dds$type,vertical=TRUE,main="SV1")
abline(h=0)
stripchart(svseq$sv[,2] ~ dds$type,vertical=TRUE,main="SV2")
abline(h=0)

ddssva <- dds
ddssva$SV1 <- svseq$sv[,1]
ddssva$SV2 <- svseq$sv[,2]
design(ddssva) <- ~ SV1 + SV2 + condition

ddssva <- DESeq(ddssva, parallel=TRUE)

ddssva$condition <- relevel(ddssva$condition, ref=ref_condition)
ressva <- results(ddssva, parallel=TRUE)
head(ressva)
# Summary of DEGs
summary(ressva)

#stabilizing the variance (takes long...)
rldssva <- rlog(ddssva)

### Removing Unknowk Varibles using RUVSeq
library(RUVSeq)
# Spikes to use as NORMALIZATOR
# Coming from literature as markers of ACTIVATION of HSC
# "LRAT","LHX2","HAND2","PDGFRB","VIM", "GFAP", "DES", "PPARG","TGFBR1","ACTA2","BAMBI"
marker_genes=c("LRAT","LHX2","HAND2","PDGFRB","VIM", "GFAP", "DES", "PPARG","TGFBR1","ACTA2","BAMBI")
# Transform to ENSID
ensgene_list = sapply(marker_genes, function(x){
  unique(genes_with_information[which(toupper(genes_with_information$mgi_symbol)==x),]$ensembl_gene_id)
})
spikes = as.vector(ensgene_list)
# Levels for condition
x = colData(dds)@listData$condition
# Extracting counts as "set" for RUV seq
set <- newSeqExpressionSet(as.matrix(raw.counts),
                           phenoData = data.frame(x, row.names=colnames(raw.counts)))

set_RUVg <- RUVg(x = set, spikes, k = 1)

plotRLE(set_RUVg, outline=FALSE, ylim=c(-4, 4))
plotPCA(set_RUVg, cex=1.2)
normalized.counts.RUV = set_RUVg@assayData$normalizedCounts
condition <- factor(rep(c("HET","WT"), each=4 ))

ddsRUV = DESeqDataSetFromMatrix(normalized.counts.RUV, DataFrame(condition), ~ condition)
dds <- DESeq(ddsRUV, parallel=TRUE)
plotDispEsts(dds)
dds<-dds[rowSums(counts(dds)) > 1,]
# Relevel the condition
dds$condition <- relevel(dds$condition, ref=ref_condition)
res <- results(dds, parallel=TRUE)
head(res)
# Summary of DEGs
summary(res)

# Stabilizing the variance (takes long...)
rld <- rlog(dds)

data_PCA <- plotPCA(rld, intgroup = c("condition","MCellCount"), returnData=TRUE)
percentVar <- round(100 * attr(data_PCA, "percentVar"))
library("ggplot2")
qplot(PC1, PC2, color=MCellCount, data=data_PCA, label=rownames(data_PCA)) +
  scale_colour_gradient(low="black", high="red") +
  geom_point(aes (size=MCellCount*5) )+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + geom_text(vjust = 1.7, nudge_x = 0.05, check_overlap = T)

library(corrplot)
col<- colorRampPalette(c("white", "blue"))(5)
corrplot(cor(normalized.counts), order = "hclust",
         hclust.method = "centroid",
         cl.lim=c(0,1),
         col= col,
         addrect = 2
         )

