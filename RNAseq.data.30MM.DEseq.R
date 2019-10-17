setwd("E:/DATI/RNAseq_30MM")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2", version = "3.8")

library(DESeq2)
library(dplyr)
library(tidyr)

RNASeq30MM.DESeq_dds <- readRDS("E:/DATI/RNAseq_30MM/RNASeq30MM.DESeq_dds.rds")


#class: DESeqDataSet 
#dim: 6 30 
#metadata(0):
#  assays(5): counts mu cooks replaceCounts replaceCooks
#rownames(6): ENSG00000223972.5 ENSG00000227232.5 ... ENSG00000238009.6 ENSG00000239945.1
#rowData names(28): baseMean baseVar ... maxCooks replace
#colnames(30): Sample_MM.021 Sample_MM.025 ... Sample_MM.406 Sample_MM.431
#colData names(44): SampleName TC ... sizeFactor replaceable

head(counts(RNASeq30MM.DESeq_dds))
dim(counts(RNASeq30MM.DESeq_dds))
#40257 features x 30 samples

#non-negative integer values in the "counts" matrix

write.table(counts(RNASeq30MM.DESeq_dds), file = "RNAseq.30MM.counts.matrix.txt", sep="\t")


colData(RNASeq30MM.DESeq_dds)
#paz.info
#DataFrame with 30 rows and 44 columns

write.table(colData(RNASeq30MM.DESeq_dds), file = "Paz.info.30MM.RNAseq.txt", sep="\t")

RNASeq30MM.DESeq_dds@assays


resultsNames(RNASeq30MM.DESeq_dds)

RNASeq30MM.DESeq_dds@design

results.trx = results(RNASeq30MM.DESeq_dds, contrast = c("Translocation", "1", "0"))
class(results.trx)


#We can order our results table by the smallest p value:
resOrdered <- results.trx[order(results.trx$pvalue),]

head(resOrdered)

#summary of the results
summary(results.trx)

#out of 40124 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 3666, 9.1%
#LFC < 0 (down)     : 2350, 5.9%
#outliers [1]       : 0, 0%
#low counts [2]     : 10244, 26%

# how many genes are significant under a FDR cut-off of 0.1
sum(results.trx$padj < 0.1, na.rm=TRUE)
#6016

#to obtain the results with FDR< 0.05
res.trx.05 <- results(RNASeq30MM.DESeq_dds, alpha=0.05)
summary(res.trx.05)

#out of 40124 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 2561, 6.4%
#LFC < 0 (down)     : 1588, 4%
#outliers [1]       : 0, 0%
#low counts [2]     : 12578, 31%


sum(res.trx.05$padj < 0.05, na.rm=TRUE)
#4149


plotMA(results.trx, ylim=c(-2,2))

#to get information about the columns of the results object
mcols(results.trx)$description


write.csv(as.data.frame(resOrdered), file="Trx_posvsneg_results_30MM.csv")




counts.30MM <-counts(RNASeq30MM.DESeq_dds, normalized=FALSE)
head(counts.30MM)

colData.30MM <- as.data.frame(colData(RNASeq30MM.DESeq_dds))
class(counts.30MM)

colData.DIS3 = colData.30MM %>% select("SampleName", "DIS3")
class(colData.DIS3)
rownames(colData.DIS3)

foo=RNASeq30MM.DESeq_dds[,-c(1:2, 4:5, 7, 9, 12)]


colData.DIS3.final = colData.DIS3[complete.cases(colData.DIS3[,2]),] 
row.names(colData.DIS3.final)
class(colData.DIS3.final)

colData.DIS3.final$condition=factor(colData.DIS3.final$DIS3)


#delete the sample columns with NA info about DIS3
counts.DIS3.MM <-counts.30MM[,-c(1:2, 4:5, 7, 9, 12)]
class(counts.DIS3.MM)
head(counts.DIS3.MM)
dim(counts.DIS3.MM)
counts.DIS3.MM


all(rownames(colData.DIS3.final) == colnames(counts.DIS3.MM))
# the analysis was run on 23 out of 30 MM samples (7 NA)
dds.DIS3 <-DESeqDataSetFromMatrix(counts.DIS3.MM, colData = colData.DIS3.final, design = ~ condition)

dds.DIS3$condition

dds.DIS3$condition <- factor(dds.DIS3$condition, levels = c("0","1"))

#running DESeq on dds.DIS3 object

dds <- DESeq(dds.DIS3)
dds <- DESeq(foo)

prova = as.matrix(read.table(file = "RNAseq.30MM.counts.matrix.txt", header = TRUE, sep = "\t"))
head(prova)
prova.paz.info=as.data.frame(read.table(file="Paz.info.30MM.RNAseq.txt", header = TRUE, sep = "\t"))
head(prova.paz.info)


colData = prova.paz.info[complete.cases(prova.paz.info$DIS3),] 
row.names(colData)
class(colData)







results.DIS3 = results(dds, contrast = c("condition", "1", "0"))

resultsNames(dds)
#[1] "Intercept"        "condition_1_vs_0"

summary(results.DIS3)
#out of 39973 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 334, 0.84%
#LFC < 0 (down)     : 81, 0.2%
#outliers [1]       : 1172, 2.9%
#low counts [2]     : 13106, 33%

sum(results.DIS3$padj < 0.1, na.rm=TRUE)
#415

res.DIS3.Ordered <- results.DIS3[order(results.DIS3$pvalue),]

write.csv(as.data.frame(res.DIS3.Ordered), file="DIS3mut_posvsneg_results_23MM.csv")


lncRNAs.DIS3corr.30MM <-read.table("RESULTS_Vanessa_Maura_lncRNAs/337lncRNA_DIS3corr_30MMRNAseq.txt", sep = "\t", header = TRUE)
dim(lncRNAs.DIS3corr.30MM)
head(lncRNAs.DIS3corr.30MM)
colnames(lncRNAs.DIS3corr.30MM)
rownames(lncRNAs.DIS3corr.30MM) <-lncRNAs.DIS3corr.30MM$probe.set
lncRNAs.DIS3corr.30MM$probe.set <-NULL


#row variance calculation across samples
#to select the first 5% (17/337) lncRNAs with the highest variance across samples
myvars.tot <- apply(lncRNAs.DIS3corr.30MM, 1, var,na.rm=TRUE) 
myvars.tot <- sort(myvars.tot,decreasing=TRUE) 

all.lncRNAs.var <-as.data.frame(myvars.tot)
head(all.lncRNAs.var)
write.table(all.lncRNAs.var, "all.lncRNAs.variance.30MM.txt", sep = "\t")
myvars <- myvars.tot[1:17]
head(myvars)
class(myvars)
str(myvars)
most.variable.5perc.lncRNAs <- as.data.frame(myvars)

lncRNAs.VAR.FILTERED <- log2(lncRNAs.DIS3corr.30MM[names(myvars),]) 
dim(lncRNAs.VAR.FILTERED) 
class(lncRNAs.VAR.FILTERED)
sample.names <- colnames(lncRNAs.DIS3corr.30MM)


library(gplots)
###HIERARCHICAL CLUSTERING BASED ON MOST VARIABLE 5% lncRNAs (based on row-variance)
scaled.lncRNAs.var.filtered <-scale(lncRNAs.VAR.FILTERED)
data.var.dist=dist(t(scaled.lncRNAs.var.filtered))
par(mfrow = c(3,1))
plot(hclust(data.var.dist),labels=sample.names, main="Complete Linkage on 5% lncRNAs with the highest variance",xlab="",sub="",ylab="")
plot(hclust(data.var.dist, method = "average"), labels=sample.names, main="Average Linkage on 5% lncRNAs with the highest variance", xlab="", ylab="")
plot(hclust(data.var.dist, method = "single"), labels=sample.names, main="Single Linkage on 5% lncRNAs with the highest variance", xlab="", ylab="")


#heatmap of 5% lncRNAs with the highest variance across 30 MM samples

par(oma=c(10,4,4,2))
my.colors = as.character(c("green", "blue","blue","blue","green","orange","green", "blue", "grey", "yellow", "yellow","yellow","grey", "yellow", "yellow","yellow", "orange", "orange", "grey","grey","blue","blue","blue","orange","grey","green","green","blue","green","green"))
ncol(lncRNAs.VAR.FILTERED)

heatmap.2(as.matrix(lncRNAs.VAR.FILTERED), Rowv = TRUE, Colv= TRUE, scale = "row", na.rm=TRUE, col=bluered, trace="none", dendrogram = "both", cexCol=1, cexRow = 0.8, ColSideColors=my.colors)

# gene lists' annotations using BioMart
library(biomaRt)
gene.list.trx = read.table("RESULTS_Vanessa_Maura_lncRNAs/52lncRNA_DIS3corr_trx.pos_vs_neg_diff.expressed_fdr1perc.txt", sep ="\t")
head(gene.list.trx)
ensglist=as.vector(gene.list.trx$V1)
length(ensglist)

ensembl = useMart( "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
geneannotation.52lncRNAs.trx <- getBM( attributes = c("ensembl_gene_id", "entrezgene", "external_gene_name", "gene_biotype"),
                         filters = "ensembl_gene_id",
                         values = ensglist,
                         mart = ensembl )

class(geneannotation.52lncRNAs.trx)
write.table(geneannotation.52lncRNAs.trx, "RESULTS_Vanessa_Maura_lncRNAs/annotation_52lncRNAs_trx_1perc.txt", sep ="\t")


gene.list.DIS3 = read.table("RESULTS_Vanessa_Maura_lncRNAs/22lncRNA_DIS3corr_DIS3mut_pos_vs_neg_diff.expressed_fdr1perc.txt", sep ="\t")
head(gene.list.DIS3)
ensglist.dis3=as.vector(gene.list.DIS3$V1)
length(ensglist.dis3)

ensembl = useMart( "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
geneannotation.22lncRNAs.dis3 <- getBM( attributes = c("ensembl_gene_id", "entrezgene", "external_gene_name", "gene_biotype"),
                                       filters = "ensembl_gene_id",
                                       values = ensglist.dis3,
                                       mart = ensembl )

listAttributes(ensembl)

class(geneannotation.22lncRNAs.dis3)
write.table(geneannotation.22lncRNAs.dis3, "RESULTS_Vanessa_Maura_lncRNAs/annotation_22lncRNAs_dis3_1perc.txt", sep ="\t")

#annotation of gene list from Gene 2.0 array analysis on DIS3MUT MM

gene.list.DIS3.Gene2 = read.table("RESULTS_Vanessa_Maura_lncRNAs/84lncRNAs_q-value0_MM_DIS3mut_20perc_vs_WT_ENSG.txt", sep ="\t")
head(gene.list.DIS3.Gene2)
ensglist.dis3.Gene2=as.vector(gene.list.DIS3.Gene2$V1)
length(ensglist.dis3.Gene2)

ensembl = useMart( "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
listAttributes(ensembl)
geneannotation.84lncRNAs.dis3.Gene2 <- getBM( attributes = c("ensembl_gene_id", "entrezgene_id","external_gene_name", "gene_biotype", "chromosome_name", "start_position","end_position","strand"),
                                        filters = "ensembl_gene_id",
                                        values = ensglist.dis3.Gene2,
                                        mart = ensembl )

write.table(geneannotation.84lncRNAs.dis3.Gene2, "RESULTS_Vanessa_Maura_lncRNAs/annotation_84lncRNAs_DIS3_qvalue0_Gene2.txt", sep ="\t")


###RUN DIFFERENTIAL EXPRESSION ANALYSIS ON DIS3mut>20% vs DIS3 WT MM
#to delete samples MM.055 and MM.310 by using their position in the list
length(colData.30MM@listData$DIS3)
idx=c(8,24)

DIS3.colData = colData.30MM@listData$DIS3[-idx]
length(DIS3.colData)
DIS3.colData


#to trasform into factor and delete NA info
condition <- as.factor(na.omit(DIS3.colData))
class(condition)


#delete the sample columns with NA info about DIS3 and with DIS3mut<20%
counts.DIS3.20perc.MM <-as.matrix(counts.30MM[,-c(1:2, 4:5, 7:8, 9, 12, 24)])
dim(counts.DIS3.20perc.MM)
all(rownames(condition) == colnames(counts.DIS3.20perc.MM))
# the analysis was run on 21 out of 30 MM samples (7 NA + 2 dis3mut<20% not included)
dds.DIS3.20perc <-DESeqDataSetFromMatrix(counts.DIS3.20perc.MM, DataFrame(condition), design = ~ condition)


#running DESeq on dds.DIS3 object

dds <- DESeq(dds.DIS3.20perc)
dds@design

results.DIS3. = results(dds, contrast = c("condition", "1", "0"))

resultsNames(dds)
#[1] "Intercept"        "condition_1_vs_0"

summary(results.DIS3)
#out of 39973 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 334, 0.84%
#LFC < 0 (down)     : 81, 0.2%
#outliers [1]       : 1172, 2.9%
#low counts [2]     : 13106, 33%

sum(results.DIS3$padj < 0.1, na.rm=TRUE)
#415

res.DIS3.Ordered <- results.DIS3[order(results.DIS3$pvalue),]

write.csv(as.data.frame(res.DIS3.Ordered), file="DIS3mut_posvsneg_results_23MM.csv")
