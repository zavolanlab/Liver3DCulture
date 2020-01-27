#!/usr/bin/env Rscript

### Created: Jan 27, 2020
### Author: Anastasiya Boersch
### Company: Zavolan Group, Biozentrum, University of Basel

# Load required libraries and set system parameters
options(java.parameters = "-Xmx10000m")
library(gplots)
library(broom)
library(ggplot2)
library(plyr)
library(grid)
library(dplyr)
library(matrixStats)
library(ggdendro)
require(gridExtra)

##### For making GO analysis with DAVID
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("RDAVIDWebService")

library("RDAVIDWebService")
splitfun <- function(a) substr(a,12,1000000L)
#####

# Load gene annotation for human
ens2name_human_rel_96 <- read.table("GeneIdTOGeneName.txt",header=F, sep="\t")
colnames(ens2name_human_rel_96) <- c("ensembl_id","gene_symbol")

# Load pre-processed data
tpms_gene <- read.table("GSE142206_tpms.txt",header=T)
rownames(tpms_gene) <- tpms_gene$ensembl_id
tpms_gene <- tpms_gene[,-c(1,2)]
tpms_gene <- tpms_gene[,c(16,17,18,1:15)]

# Annotation of replicates
des=c("Control","Collagen (0.5 mg/ml)","Collagen (1.0 mg/ml)","Matrigel (1.5 mg/ml)","Matrigel (3.0 mg/ml)","Matrigel (6.0 mg/ml)")
n_rep <- c(3,3,3,3,3,3)

# Plot coordinates of nth (n_comp) component
plot_pca <- function(svd.mat,n,des,n_rep,n_comp){
  n_cond=length(n_rep)
  ev <- svd.mat$v
  d_2 <- (svd.mat$d)^2
  par(mar=c(8.1, 4.1, 4.1, 4.1), xpd=TRUE)
  df_ggplot <- data.frame(des=rep(des,n_rep),coord=ev[((n_comp-1)*n+1):(n_comp*n)])
  df_ggplot$des <- factor(df_ggplot$des, levels = des)
  g <- (ggplot(df_ggplot, aes(x=des, y=coord)) + 
          geom_boxplot(outlier.shape = NA) +
          geom_jitter(shape=19, position=position_jitter(0.2),size=2.5) +
          theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                             text = element_text(family="Helvetica",size=8,color="black"),axis.title.x=element_blank(),
                             axis.text.x = element_text(family="Helvetica",size=8,angle = 45,hjust = 1),
                             axis.title.y = element_text(family="Helvetica",color="black", size=9)) +
          ylab(paste("PC",toString(n_comp)," (explained variance: ",round(d_2[n_comp]/sum(d_2)*100,3),"%)",sep="")))
  return(g)
}

# Get genes contributing to a principle component
pattern <- function(svd_mat,gene_list,n_rep,n_comp) {
  ev <- svd_mat$v
  # Calculate projections
  proj <- svd_mat$u%*%diag(svd_mat$d)
  rownames(proj) <- gene_list
  # Calculate correlations
  crl <- proj/sqrt(rowSums(proj^2))
  # Scaling projections
  sc.proj <- proj[,n_comp]/max(abs(proj[,n_comp]))
  # Create a data frame with projections, correlations and z scores for projections
  df <- data.frame(proj=sc.proj,cor=crl[,n_comp],proj.z = (sc.proj-mean(sc.proj))/sd(sc.proj))
  rownames(df) <- rownames(proj)
  df$ensembl_id <- rownames(df)
  return(df)
}

# Perform principal component analysis for all conditions together
# Remove lowly expressed transcripts and perform SVD
keep <- rowSums(tpms_gene>1)>=min(n_rep)
tpms_gene_k <- tpms_gene[keep,]
tpms_gene.log <- log2(tpms_gene_k+1)

tpms_gene_s1 <- apply(tpms_gene.log,2,scale,scale=F,center=T)
tpms_gene_s2 <- apply(tpms_gene_s1,1,scale,scale=F,center=T)

svd.mat <- svd(t(tpms_gene_s2))
ev <- svd.mat$v
d <- svd.mat$d

pc1 <- ev[,1] # Increasing profile
pc1_var <- round(d[1]^2/sum(d^2)*100,digits=2)
pc2_var <- round(d[2]^2/sum(d^2)*100,digits=2)

plot_df <- data.frame(PC1=pc1,PC2=ev[,2])
plot_df$cond <- factor(rep(des,n_rep),levels = des)
p1 <- ggplot(plot_df, aes(x=PC1, y=PC2, color=cond))+
  geom_point()+
  scale_color_manual(values=c("cyan", "blue", "green","red","black","orange"))+
  theme_bw()+
  coord_fixed()+
  xlab(paste0("PC1, ",pc1_var,"%"))+
  ylab(paste0("PC2, ",pc2_var,"%"))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.text.x = element_text(family="Helvetica",color="black", size=8),
        axis.text.y = element_text(family="Helvetica",color="black", size=8),
        axis.title.x = element_text(family="Helvetica",color="black", size=9),
        axis.title.y = element_text(family="Helvetica",color="black", size=9))
p1

p1=plot_pca(svd.mat,sum(n_rep),des,n_rep,1)
p2=plot_pca(svd.mat,sum(n_rep),des,n_rep,2)

grid.newpage()
grid.draw(cbind(ggplotGrob(p1), ggplotGrob(p2),size = "last"))

# Align gene expression to PCs
pattern_genes_pc1 <- pattern(svd.mat,rownames(tpms_gene_k),n_rep,1)
pattern_genes_pc1 <- merge(ens2name_human_rel_96,pattern_genes_pc1,by="ensembl_id")
ord <- order(pattern_genes_pc1$proj)
pattern_genes_pc1 <- pattern_genes_pc1[ord,]

pattern_genes_pc2 <- pattern(svd.mat,rownames(tpms_gene_k),n_rep,2)
pattern_genes_pc2 <- merge(ens2name_human_rel_96,pattern_genes_pc2,by="ensembl_id")
ord <- order(pattern_genes_pc2$proj)
pattern_genes_pc2 <- pattern_genes_pc2[ord,]

# Subset genes
zscore_thresh <- 1.96
cor_thresh <- 0.6
pc1_genes_higher_cg <- pattern_genes_pc1[which(pattern_genes_pc1$proj.z>=zscore_thresh & pattern_genes_pc1$cor>=cor_thresh),]
pc1_genes_lower_cg <- pattern_genes_pc1[which(pattern_genes_pc1$proj.z<=-zscore_thresh & pattern_genes_pc1$cor<=-cor_thresh),]

pc2_genes_higher_mg <- pattern_genes_pc2[which(pattern_genes_pc2$proj.z>=zscore_thresh & pattern_genes_pc2$cor>=cor_thresh),]
pc2_genes_lower_mg <- pattern_genes_pc2[which(pattern_genes_pc2$proj.z<=-zscore_thresh & pattern_genes_pc2$cor<=-cor_thresh),]

# Make a heatmap for genes aligned with PC1 and PC2
df <- tpms_gene.log[c(as.character(pc1_genes_higher_cg$ensembl_id),
                      as.character(pc1_genes_lower_cg$ensembl_id),
                      as.character(pc2_genes_higher_mg$ensembl_id),
                      as.character(pc2_genes_lower_mg$ensembl_id)),]
df <- (as.matrix(df)-rowMeans(as.matrix(df)))/rowSds(as.matrix(df))
df[df > 1] <- 1
df[df < -1] <- -1
df <- df[nrow(df):1,]
heatmap_ggplot <- expand.grid(gene_ids = rownames(df),replicates = colnames(df))
heatmap_ggplot$zsc <- as.vector(as.matrix(df))

n_tot <- nrow(pc1_genes_higher_cg)+nrow(pc1_genes_lower_cg)+nrow(pc2_genes_higher_mg)+nrow(pc2_genes_lower_mg)
y1 <- nrow(pc2_genes_lower_mg)
y2 <- y1+nrow(pc2_genes_higher_mg)
y3 <- y2+nrow(pc1_genes_lower_cg)
my.lines<-data.frame(x=c(3,6,9,12,15,0,0,0)+0.5,y=c(0,0,0,0,0,y1,y2,y3)+0.5,xend=c(3,6,9,12,15,18,18,18)+0.5,yend=c(n_tot,n_tot,n_tot,n_tot,n_tot,y1,y2,y3)+0.5)

p1 <- (ggplot(data=heatmap_ggplot, aes(x = replicates, y = gene_ids)) 
  + geom_tile(aes(fill = zsc))
  # indicate paramters of the colobar of the legend
  + scale_fill_gradient2(low = "blue", high = "red", mid = "white",midpoint = 0, space = "Lab",name="mRNA level (z-score)")
  # put the legend on the top 
  + theme(legend.position="top",
          axis.ticks.y = element_blank(),
          legend.text=element_text(family="Helvetica",color="black", size=8),
          axis.title.x=element_blank(),axis.title.y=element_blank(),
          axis.text.x = element_text(family="Helvetica",size=8,angle=45,hjust = 1,color = "black"),
          axis.text.y = element_text(family="Helvetica",color="black", size=8))
  # put gene labels on the right and reverse the order
  + scale_y_discrete(position="right",labels=c())
  + scale_x_discrete(labels=c("Control, rep1","Control, rep2","Control, rep3",
                              "Collagen (0.5 mg/ml), rep1","Collagen (0.5 mg/ml), rep2","Collagen (0.5 mg/ml), rep3",
                              "Collagen (1.0 mg/ml), rep1","Collagen (1.0 mg/ml), rep2","Collagen (1.0 mg/ml), rep3",
                              "Matrigel (1.5 mg/ml), rep1","Matrigel (1.5 mg/ml), rep2","Matrigel (1.5 mg/ml), rep3",
                              "Matrigel (3.0 mg/ml), rep1","Matrigel (3.0 mg/ml), rep2","Matrigel (3.0 mg/ml), rep3",
                              "Matrigel (6.0 mg/ml), rep1","Matrigel (6.0 mg/ml), rep2","Matrigel (6.0 mg/ml), rep3"))
  + geom_segment(data=my.lines, aes(x,y,xend=xend, yend=yend), size=0.5, inherit.aes=F)
  + theme(axis.title.x=element_blank(),axis.title.y=element_blank()))

print(p1)

# Perform DAVID analysis of genes aligned with PC1
david<-DAVIDWebService(email="your_email", url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
# Problems with timing out while loading the background https://support.bioconductor.org/p/65696/
setTimeOut(david, 100000)
setAnnotationCategories(david, c("GOTERM_BP_DIRECT","GOTERM_MF_DIRECT", "GOTERM_CC_DIRECT"))

# Load background
BG <- addList(david,rownames(tpms_gene_k),idType = "ENSEMBL_GENE_ID",listName = "expr_genes",listType = "Background")

# Annotation of genes higher in CG than in Control
FG <- addList(david,as.character(pc1_genes_higher_cg$ensembl_id),idType = "ENSEMBL_GENE_ID",listName = "pc1_genes_higher_cg",listType = "Gene")
setCurrentBackgroundPosition(david, 2)
setCurrentGeneListPosition(david,1)
david
david_pc1_genes_higher_cg <- data.frame(getFunctionalAnnotationChart(david))

# Annotation of genes lower in CG than in Control
FG <- addList(david,as.character(pc1_genes_lower_cg$ensembl_id),idType = "ENSEMBL_GENE_ID",listName = "pc1_genes_lower_cg",listType = "Gene")
setCurrentBackgroundPosition(david, 2)
setCurrentGeneListPosition(david,2)
david
david_pc1_genes_lower_cg <- data.frame(getFunctionalAnnotationChart(david))

# Make barplots for top 20 enriched terms
# Barplot for top 20 GO terms for genes positively correlated with PC1
ax_lim <- min(c(log10(david_pc1_genes_lower_cg[20:1,"PValue"]),log10(david_pc1_genes_higher_cg[20:1,"PValue"])))
go_names <- as.character(david_pc1_genes_higher_cg[20:1,"Term"])
go_names <- apply(as.matrix(go_names),2,splitfun)
df <- data.frame(name=go_names,log10_pv=log10(david_pc1_genes_higher_cg[20:1,"PValue"]))
df$name <- factor(df$name, levels = df$name)
p1 <- ggplot(df, aes(name, log10_pv)) +
  geom_col() +
  coord_flip() +
  scale_y_reverse(limits=c(0,ax_lim))+
  scale_x_discrete(position="top") +
  theme_bw() +
  theme(axis.text.x = element_text(family="Helvetica",colour = "black", size=8),axis.text.y = element_text(family="Helvetica",colour = "black", size=8),
        axis.title.y=element_blank(),axis.title.x=element_text(family="Helvetica",colour = "black", size=9),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_abline(intercept=2, slope=0, linetype=2, color="red", size=0.5)+
  ylab("log10(p-value)")

# Barplot for top 20 GO terms for genes negatively correlated with PC1
go_names <- as.character(david_pc1_genes_lower_cg[20:1,"Term"])
go_names <- apply(as.matrix(go_names),2,splitfun)
df <- data.frame(name=go_names,log10_pv=log10(david_pc1_genes_lower_cg[20:1,"PValue"]))
df$name <- factor(df$name, levels = df$name)
p2 <- ggplot(df, aes(name, log10_pv)) +
  geom_col() +
  coord_flip() +
  ylim(ax_lim,0) +
  theme_bw() +
  theme(axis.text.x = element_text(family="Helvetica",colour = "black", size=8),axis.text.y = element_text(family="Helvetica",colour = "black", size=8),
        axis.title.y=element_blank(),axis.title.x=element_text(family="Helvetica",colour = "black", size=9),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_abline(intercept=-2, slope=0, linetype=2, color="red", size=0.5)+
  ylab("log10(p-value)") 

grid.newpage()
grid.draw(cbind(ggplotGrob(p2), ggplotGrob(p1),size = "last"))

# Get the list of genes contributing to extracellular region/extreacellular space; genes were used for the STRING annotation
tmp1 <- as.character(unlist(strsplit(david_pc1_genes_higher_cg[1,6],", ")))
tmp2 <- as.character(unlist(strsplit(david_pc1_genes_lower_cg[1,6],", ")))
View(data.frame(unique(c(tmp1,tmp2))))

# Perform DAVID analysis of genes aligned with PC2

# Annotation of genes higher in CG than in Control
FG <- addList(david,as.character(pc2_genes_higher_mg$ensembl_id),idType = "ENSEMBL_GENE_ID",listName = "pc2_genes_higher_mg",listType = "Gene")
setCurrentBackgroundPosition(david, 2)
setCurrentGeneListPosition(david,3)
david
david_pc2_genes_higher_mg <- data.frame(getFunctionalAnnotationChart(david))

# Annotation of genes lower in CG than in Control
FG <- addList(david,as.character(pc2_genes_lower_mg$ensembl_id),idType = "ENSEMBL_GENE_ID",listName = "pc2_genes_lower_mg",listType = "Gene")
setCurrentBackgroundPosition(david, 2)
setCurrentGeneListPosition(david,4)
david
david_pc2_genes_lower_mg <- data.frame(getFunctionalAnnotationChart(david))

# Make barplots for top 20 enriched terms
# Barplot for top 20 GO terms for genes positively correlated with PC2
ax_lim <- min(c(log10(david_pc2_genes_lower_mg[20:1,"PValue"]),log10(david_pc2_genes_higher_mg[20:1,"PValue"])))
go_names <- as.character(david_pc2_genes_higher_mg[20:1,"Term"])
go_names <- apply(as.matrix(go_names),2,splitfun)
df <- data.frame(name=go_names,log10_pv=log10(david_pc2_genes_higher_mg[20:1,"PValue"]))
df$name <- factor(df$name, levels = df$name)
p1 <- ggplot(df, aes(name, log10_pv)) +
  geom_col() +
  coord_flip() +
  scale_y_reverse(limits=c(0,ax_lim))+
  scale_x_discrete(position="top") +
  theme_bw() +
  theme(axis.text.x = element_text(family="Helvetica",colour = "black", size=8),axis.text.y = element_text(family="Helvetica",colour = "black", size=8),
        axis.title.y=element_blank(),axis.title.x=element_text(family="Helvetica",colour = "black", size=9),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_abline(intercept=2, slope=0, linetype=2, color="red", size=0.5)+
  ylab("log10(p-value)")

# Barplot for top 20 GO terms for genes negatively correlated with PC2
go_names <- as.character(david_pc2_genes_lower_mg[20:1,"Term"])
go_names <- apply(as.matrix(go_names),2,splitfun)
df <- data.frame(name=go_names,log10_pv=log10(david_pc2_genes_lower_mg[20:1,"PValue"]))
df$name <- factor(df$name, levels = df$name)
p2 <- ggplot(df, aes(name, log10_pv)) +
  geom_col() +
  coord_flip() +
  ylim(ax_lim,0) +
  theme_bw() +
  theme(axis.text.x = element_text(family="Helvetica",colour = "black", size=8),axis.text.y = element_text(family="Helvetica",colour = "black", size=8),
        axis.title.y=element_blank(),axis.title.x=element_text(family="Helvetica",colour = "black", size=9),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_abline(intercept=-2, slope=0, linetype=2, color="red", size=0.5)+
  ylab("log10(p-value)") 

grid.newpage()
grid.draw(cbind(ggplotGrob(p2), ggplotGrob(p1),size = "last"))

# Get the list of genes contributing to extracellular region/extreacellular space; genes were used for the STRING annotation
tmp1 <- as.character(unlist(strsplit(david_pc2_genes_higher_mg[1,6],", ")))
tmp2 <- as.character(unlist(strsplit(david_pc2_genes_lower_mg[1,6],", ")))
View(data.frame(unique(c(tmp1,tmp2))))





























