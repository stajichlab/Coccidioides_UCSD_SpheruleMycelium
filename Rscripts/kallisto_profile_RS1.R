library(DESeq2)
library(tximport)
library(dplyr)
library(ggplot2)
library(magrittr)
library(Biobase)
library(pheatmap)
library(RColorBrewer)
library(IHW)
library(apeglm)
samples <- read.csv("samples.csv",header=TRUE)
exprnames <- do.call(paste,c(samples[c("Condition","Replicate")],sep="."))
exprnames <- sub(".([123])$",".r\\1",exprnames,perl=TRUE)
files <- file.path("results","RS1",exprnames,"abundance.h5")
txi.kallisto <- tximport(files, type = "kallisto", txOut = TRUE)
head(txi.kallisto$counts)
colnames(txi.kallisto$counts) <- exprnames
colnames(txi.kallisto$abundance) <- exprnames
write.csv(txi.kallisto$abundance,"reports/RS1_kallisto.TPM.csv")
write.csv(txi.kallisto$counts,"reports/RS1_kallisto.counts.csv")

# DEseq2 analyses
treatment = factor (samples$Condition)

sampleTable <- data.frame(condition = treatment)
rownames(sampleTable) = exprnames

dds <- DESeqDataSetFromTximport(txi.kallisto,sampleTable,design=~ condition)

#nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ]
#nrow(dds)

dds <- estimateSizeFactors(dds)
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
#head(assay(vsd), 3)

dds <- DESeq(dds,fitType='mean')
resultsNames(dds)
#dds <- lfcShrink(dds,coef="condition_Spherule48H_vs_Mycelia",type="apeglm")
#res <- results(dds, contrast=c("condition","Mycelia","Spherule48H"))
#res <- results(dds, name="condition_Spherule48h_vs_Mycelia")
res <- results(dds,filterFun=ihw)
summary(res)
res <- subset(res,res$padj<0.05)
res <- res[order(res$pvalue ),]
summary(res)


df <- bind_rows(
  as_tibble(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
         mutate(transformation = "log2(x + 1)"),
  as_tibble(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"),
  as_tibble(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"))

colnames(df)[1:2] <- c("x", "y")

pdf("plots/RS1_RNASeq_kallisto.pdf")
ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)

select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:75]
df <- as.data.frame(colData(dds)[,c("condition")])
rownames(df) = exprnames
colnames(df) = c("Condition")

pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df,main="VSD ordered")

topVar <- head(order(rowVars(assay(vsd)),
                                decreasing=TRUE),60)
mat  <- assay(vsd)[ topVar, ]
mat  <- mat - rowMeans(mat)
pheatmap(mat, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df,main="VSD - most different")

pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df,main="RLD")

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=condition)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    coord_fixed()

## res <- results(dds, contrast=c("condition","Mycelia","Spherule48H"))
## summary(res)
## plotMA(res,main="Mycelium vs Spherule48H")
## res <- subset(res,res$padj<0.05)
## res <- res[order(res$pvalue ),]
## summary(res)
## write.csv(res,"reports/kallisto.DESeq_Mycelium_Spherule48H.csv")


## res <- results(dds, contrast=c("condition","Mycelia","Spherule8D"))
## summary(res)
## plotMA(res,main="Mycelium vs Spherule8D")
## res <- subset(res,res$padj<0.05)
## res <- res[order(res$pvalue ),]
## summary(res)
## write.csv(res,"reports/kallisto.DESeq_Mycelium_Spherule8D.csv")

## res <- results(dds, contrast=c("condition","Spherule48H","Spherule8D"))
## summary(res)
## plotMA(res,main="Spherule48H vs Spherule8D")
## res <- subset(res,res$padj<0.05)
## res <- res[order(res$pvalue ),]
## summary(res)
## write.csv(res,"reports/kallisto.DESeq_Spherule48H_Spherule8D.csv")
