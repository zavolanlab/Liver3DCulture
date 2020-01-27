#!/usr/bin/env Rscript

### Created: Jan 27, 2020
### Author: Anastasiya Boersch
### Company: Zavolan Group, Biozentrum, University of Basel

# Load required libraries
#source("http://bioconductor.org/biocLite.R")
library(gplots)
library(edgeR)
library(broom)
library(ggplot2)
library(plyr)
library(grid)
library(dplyr)
library(matrixStats)

# Load gene annotation for human
ens2name_human_rel_96 <- read.table("GeneIdTOGeneName.txt",header=F, sep="\t")
colnames(ens2name_human_rel_96) <- c("ensembl_id","gene_symbol")

# Load pre-processed data
counts_gene_rnaseq <- read.table("GSE142206_counts.txt",header=T)
rownames(counts_gene_rnaseq) <- as.character(counts_gene_rnaseq$ensembl_id)
counts_gene_rnaseq <- counts_gene_rnaseq[,-c(1,2)]

# Annotation of replicates
des=c("CG_0.5","CG_1.0","MG_1.5","MG_3.0","MG_6.0","CON")
n_rep <- c(3,3,3,3,3,3)

# DE analysis for conditions
# Annotate conditions to replicates
group <- factor(rep(des,n_rep),levels = des)
# DE analysis is based on the analysis of raw counts
cts <- round(counts_gene_rnaseq)
y <- DGEList(counts=cts,group=group)
# Normalization
y <- calcNormFactors(y)

# Quantify counts per million
cpm_ent <- data.frame(cpm(y))

# Filtering lowly expressed genes
keep <- rowSums(cpm(y)>1)>=min(n_rep)
y <- y[keep,,keep.lib.sizes=F]

# Create a design matrix and estimate dispersion
design <- model.matrix(~group)
y <- estimateDisp(y,design,robust=TRUE)
plotBCV(y)

# Testing for DE genes
fit <- glmFit(y,design)

# Initialize matrices for storing the results of all pairwise comparisons
# Matrix with log fold changes
DE_logFC <- matrix(ncol=(length(des)/2)*(length(des)-1), nrow=nrow(y$counts))
# Matrix with p-values
DE_pval <- matrix(ncol=(length(des)/2)*(length(des)-1), nrow=nrow(y$counts))
# Matrix with FDRs
DE_FDR <- matrix(ncol=(length(des)/2)*(length(des)-1), nrow=nrow(y$counts))
# Matrix with median expression per condition
DE_expr <- matrix(ncol=(length(des)/2)*(length(des)-1), nrow=nrow(y$counts))
coln <- c()
# Perform comparisons of all conditions one by one to the 1st one
k=1
for (j in seq(2,length(des))){
  # This is the difference between this and next for loops
  DE <- glmLRT(fit,coef=j)
  DEt <- DE$table
  DE_logFC[,k]=DEt$logFC
  DE_pval[,k]=DEt$PValue
  DE_expr[,k]=DEt$logCPM
  DE_FDR[,k]=p.adjust(DEt$PValue,method="BH")
  k=k+1
  coln <- c(coln,paste(des[j],"v",des[1],sep=""))
}
for (i in seq(2,length(des)-1)){
  for (j in seq(i+1,length(des))){
    # Create a constrast for comparing two conditions between each other
    x <- rep(0,length(des))
    x[i]=-1
    x[j]=1
    # This is the difference between this and previous for loops
    DE <- glmLRT(fit,contrast=x)
    DEt <- DE$table
    DE_logFC[,k]=DEt$logFC
    DE_pval[,k]=DEt$PValue
    DE_expr[,k]=DEt$logCPM
    DE_FDR[,k]=p.adjust(DEt$PValue,method="BH")
    k=k+1
    coln <- c(coln,paste(des[j],"v",des[i],sep=""))
  }
}

DE_logFC <- data.frame(DE_logFC)
DE_pval <- data.frame(DE_pval)
DE_FDR <- data.frame(DE_FDR)
DE_expr <- data.frame(DE_expr)
rownames(DE_logFC) <- rownames(y$counts)
rownames(DE_pval) <- rownames(y$counts)
rownames(DE_FDR) <- rownames(y$counts)
rownames(DE_expr) <- rownames(y$counts)
colnames(DE_logFC) <- coln
colnames(DE_pval) <- coln
colnames(DE_FDR) <- coln
colnames(DE_expr) <- coln

# Create a directory for saving results of DE analysis
system("mkdir -p DE_analysis")

# Save the resulting table with the DE analysis for sharing
DE_logFC_save <- DE_logFC
DE_logFC_save$ensembl_id <- rownames(DE_logFC_save)
DE_logFC_save <- merge(ens2name_human_rel_96,DE_logFC_save,by="ensembl_id")
write.table(DE_logFC_save,"DE_analysis/DE_logFC.txt",quote = F,col.names = T,row.names = F,sep="\t")

DE_FDR_save <- DE_FDR
DE_FDR_save$ensembl_id <- rownames(DE_FDR_save)
DE_FDR_save <- merge(ens2name_human_rel_96,DE_FDR_save,by="ensembl_id")
write.table(DE_FDR_save,"DE_analysis/DE_FDR.txt",quote = F,col.names = T,row.names = F,sep="\t")




