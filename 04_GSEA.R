#!/usr/bin/env Rscript

### Created: Jan 27, 2020
### Author: Anastasiya Boersch
### Company: Zavolan Group, Biozentrum, University of Basel

# Load required libraries
library(grid)
library(ggplot2)

# GSEA
# Download GSEA for command line and follow instructions, see https://software.broadinstitute.org/gsea/downloads.jsp

# Function for runnung GSEA
gsea_command_line <- function(gmt_file,rnk_file,report_label,out_dir) {
  system(paste("ml purge &&","ml Java/11.0.3_7 &&",'export _JAVA_OPTIONS="-Xmx1024m" &&', "sh gsea_software/GSEA_4.0.3/gsea-cli.sh GSEAPreranked",
               "-gmx", gmt_file,
               "-rnk", rnk_file,
               "-set_max 500",
               "-set_min 15",
               "-plot_top_x 50",
               "-rpt_label", report_label,
               "-zip_report false",
               "-out", out_dir,sep=" "),wait=FALSE)
}

# Prepare tables for GSEA
# Load human KEGG pathways and genes assigned to pathways for human
path_id_name_human <- read.table("HumanPathIdName.txt",sep="\t")
path_gene_human <- read.table("HumanGenesPathways.txt",sep="\t")

# Remove pathways associated with deseases
path_id_name_human <- path_id_name_human[which(substr(path_id_name_human$path_id,10,10)!="5"),]
path_gene_human <- path_gene_human[which(substr(path_gene_human$path_id,10,10)!="5"),]
path_id_name_human <- path_id_name_human[-which(path_id_name_human$path_id  %in% c("path:hsa04930","path:hsa04940","path:hsa04950","path:hsa04932","path:hsa04931","path:hsa04933","path:hsa04934")),]
path_gene_human <- path_gene_human[-which(path_gene_human$path_id  %in% c("path:hsa04930","path:hsa04940","path:hsa04950","path:hsa04932","path:hsa04931","path:hsa04933","path:hsa04934")),]

# Create a table containing the code of the pathway and corresponding genes (gmt format required for GSEA)
mat_path <- matrix(,ncol=1500,nrow=length(path_id_name_human$path_id))
k <- 1
for (i in as.character(path_id_name_human$path_id)){
  mat_path[k,1] <- i
  genes_path <- as.character(path_gene_human[which(path_gene_human$path_id == i),"gene_symbol"])
  j=3
  for (g in genes_path){
    mat_path[k,j] <- g
    j <- j+1
  }
  k <- k+1
}
mat_path[1:(k-1),2] <- rep("na",k-1)

# Create a directory for storing GSEA
system("mkdir -p GSEA")
write.table(mat_path,"GSEA/KEGG_pathways_human.gmt",sep="\t",row.names = F,col.names = F,quote = F, na="")

# Load the table from DE analysis containing log fold changes
DE_logFC <- read.table("DE_analysis/DE_logFC.txt",sep="\t",header = T)

##### GSEA for Control vs MG 1.5
# Create a repository for storing the output
system("mkdir -p GSEA/GSEA_CONvMG_1.5")

# Prepare ranking tables for Control vs MG 1.5
DE_logFC_GSEA <- DE_logFC[,c("gene_symbol","CONvMG_1.5")]
ord <- order(DE_logFC_GSEA$CONvMG_1.5,decreasing = T)
DE_logFC_GSEA <- DE_logFC_GSEA[ord,]
write.table(DE_logFC_GSEA,"GSEA/GSEA_CONvMG_1.5/Preranked_CONvMG_1.5.rnk",sep="\t",row.names = F,col.names = F,quote = F)

# Specify the location of required repositories and files
gmt_file <- "GSEA/KEGG_pathways_human.gmt"
rnk_file <- "GSEA/GSEA_CONvMG_1.5/Preranked_CONvMG_1.5.rnk"
report_label <- "GSEA_CONvMG_1.5"
out_dir <- "GSEA/GSEA_CONvMG_1.5"
gsea_command_line(gmt_file,rnk_file,report_label,out_dir)

##### GSEA for Control vs MG 3
# Create a repository for storing the output
system("mkdir -p GSEA/GSEA_CONvMG_3.0")

# Prepare ranking tables for Control vs MG 3
DE_logFC_GSEA <- DE_logFC[,c("gene_symbol","CONvMG_3.0")]
ord <- order(DE_logFC_GSEA$CONvMG_3.0,decreasing = T)
DE_logFC_GSEA <- DE_logFC_GSEA[ord,]
write.table(DE_logFC_GSEA,"GSEA/GSEA_CONvMG_3.0/Preranked_CONvMG_3.0.rnk",sep="\t",row.names = F,col.names = F,quote = F)

# Specify the location of required repositories and files
gmt_file <- "GSEA/KEGG_pathways_human.gmt"
rnk_file <- "GSEA/GSEA_CONvMG_3.0/Preranked_CONvMG_3.0.rnk"
report_label <- "GSEA_CONvMG_3.0"
out_dir <- "GSEA/GSEA_CONvMG_3.0"
gsea_command_line(gmt_file,rnk_file,report_label,out_dir)

##### GSEA for Control vs MG 6
# Create a repository for storing the output
system("mkdir -p GSEA/GSEA_CONvMG_6.0")

# Prepare ranking tables for Control vs MG 6
DE_logFC_GSEA <- DE_logFC[,c("gene_symbol","CONvMG_6.0")]
ord <- order(DE_logFC_GSEA$CONvMG_6.0,decreasing = T)
DE_logFC_GSEA <- DE_logFC_GSEA[ord,]
write.table(DE_logFC_GSEA,"GSEA/GSEA_CONvMG_6.0/Preranked_CONvMG_6.0.rnk",sep="\t",row.names = F,col.names = F,quote = F)

# Specify the location of required repositories and files
gmt_file <- "GSEA/KEGG_pathways_human.gmt"
rnk_file <- "GSEA/GSEA_CONvMG_6.0/Preranked_CONvMG_6.0.rnk"
report_label <- "GSEA_CONvMG_6.0"
out_dir <- "GSEA/GSEA_CONvMG_6.0"
gsea_command_line(gmt_file,rnk_file,report_label,out_dir)

##### GSEA for Control vs CG 0.5
# Create a repository for storing the output
system("mkdir -p GSEA/GSEA_CONvCG_0.5")

# Prepare ranking tables for Control vs CG 0.5
DE_logFC_GSEA <- DE_logFC[,c("gene_symbol","CONvCG_0.5")]
ord <- order(DE_logFC_GSEA$CONvCG_0.5,decreasing = T)
DE_logFC_GSEA <- DE_logFC_GSEA[ord,]
write.table(DE_logFC_GSEA,"GSEA/GSEA_CONvCG_0.5/Preranked_CONvCG_0.5.rnk",sep="\t",row.names = F,col.names = F,quote = F)

# Specify the location of required repositories and files
gmt_file <- "GSEA/KEGG_pathways_human.gmt"
rnk_file <- "GSEA/GSEA_CONvCG_0.5/Preranked_CONvCG_0.5.rnk"
report_label <- "GSEA_CONvCG_0.5"
out_dir <- "GSEA/GSEA_CONvCG_0.5"
gsea_command_line(gmt_file,rnk_file,report_label,out_dir)

##### GSEA for Control vs CG 1
# Create a repository for storing the output
system("mkdir -p GSEA/GSEA_CONvCG_1.0")

# Prepare ranking tables for Control vs CG 1
DE_logFC_GSEA <- DE_logFC[,c("gene_symbol","CONvCG_1.0")]
ord <- order(DE_logFC_GSEA$CONvCG_1.0,decreasing = T)
DE_logFC_GSEA <- DE_logFC_GSEA[ord,]
write.table(DE_logFC_GSEA,"GSEA/GSEA_CONvCG_1.0/Preranked_CONvCG_1.0.rnk",sep="\t",row.names = F,col.names = F,quote = F)

# Specify the location of required repositories and files
gmt_file <- "GSEA/KEGG_pathways_human.gmt"
rnk_file <- "GSEA/GSEA_CONvCG_1.0/Preranked_CONvCG_1.0.rnk"
report_label <- "GSEA_CONvCG_1.0"
out_dir <- "GSEA/GSEA_CONvCG_1.0"
gsea_command_line(gmt_file,rnk_file,report_label,out_dir)

# Load GSEA reports !!!!!!!!! Note that names of output directories with GSEA results contain a time stamp !!!!!!!!!!!!!
# Control vs MG 1.5
gsea_con_v_mg1.5_pos <- read.table("GSEA/GSEA_CONvMG_1.5/GSEA_CONvMG_1.5.GseaPreranked.1578923228839/gsea_report_for_na_pos_1578923228839.xls",sep="\t")
gsea_con_v_mg1.5_pos <- gsea_con_v_mg1.5_pos[,-c(2,3,12)]
colnames(gsea_con_v_mg1.5_pos) <- unlist(gsea_con_v_mg1.5_pos[1,])
gsea_con_v_mg1.5_pos <- gsea_con_v_mg1.5_pos[-1,]

gsea_con_v_mg1.5_neg <- read.table("GSEA/GSEA_CONvMG_1.5/GSEA_CONvMG_1.5.GseaPreranked.1578923228839/gsea_report_for_na_neg_1578923228839.xls",sep="\t")
gsea_con_v_mg1.5_neg <- gsea_con_v_mg1.5_neg[,-c(2,3,12)]
colnames(gsea_con_v_mg1.5_neg) <- unlist(gsea_con_v_mg1.5_neg[1,])
gsea_con_v_mg1.5_neg <- gsea_con_v_mg1.5_neg[-1,]

gsea_con_v_mg1.5 <- rbind(gsea_con_v_mg1.5_pos,gsea_con_v_mg1.5_neg)

# Control vs MG 3.0
gsea_con_v_mg3.0_pos <- read.table("GSEA/GSEA_CONvMG_3.0/GSEA_CONvMG_3.0.GseaPreranked.1578923228843/gsea_report_for_na_pos_1578923228843.xls",sep="\t")
gsea_con_v_mg3.0_pos <- gsea_con_v_mg3.0_pos[,-c(2,3,12)]
colnames(gsea_con_v_mg3.0_pos) <- unlist(gsea_con_v_mg3.0_pos[1,])
gsea_con_v_mg3.0_pos <- gsea_con_v_mg3.0_pos[-1,]

gsea_con_v_mg3.0_neg <- read.table("GSEA/GSEA_CONvMG_3.0/GSEA_CONvMG_3.0.GseaPreranked.1578923228843/gsea_report_for_na_neg_1578923228843.xls",sep="\t")
gsea_con_v_mg3.0_neg <- gsea_con_v_mg3.0_neg[,-c(2,3,12)]
colnames(gsea_con_v_mg3.0_neg) <- unlist(gsea_con_v_mg3.0_neg[1,])
gsea_con_v_mg3.0_neg <- gsea_con_v_mg3.0_neg[-1,]

gsea_con_v_mg3.0 <- rbind(gsea_con_v_mg3.0_pos,gsea_con_v_mg3.0_neg)

# Control vs MG 6.0
gsea_con_v_mg6.0_pos <- read.table("GSEA/GSEA_CONvMG_6.0/GSEA_CONvMG_6.0.GseaPreranked.1578923228842/gsea_report_for_na_pos_1578923228842.xls",sep="\t")
gsea_con_v_mg6.0_pos <- gsea_con_v_mg6.0_pos[,-c(2,3,12)]
colnames(gsea_con_v_mg6.0_pos) <- unlist(gsea_con_v_mg6.0_pos[1,])
gsea_con_v_mg6.0_pos <- gsea_con_v_mg6.0_pos[-1,]

gsea_con_v_mg6.0_neg <- read.table("GSEA/GSEA_CONvMG_6.0/GSEA_CONvMG_6.0.GseaPreranked.1578923228842/gsea_report_for_na_neg_1578923228842.xls",sep="\t")
gsea_con_v_mg6.0_neg <- gsea_con_v_mg6.0_neg[,-c(2,3,12)]
colnames(gsea_con_v_mg6.0_neg) <- unlist(gsea_con_v_mg6.0_neg[1,])
gsea_con_v_mg6.0_neg <- gsea_con_v_mg6.0_neg[-1,]

gsea_con_v_mg6.0 <- rbind(gsea_con_v_mg6.0_pos,gsea_con_v_mg6.0_neg)

# Control vs CG 0.5
gsea_con_v_cg0.5_pos <- read.table("GSEA/GSEA_CONvCG_0.5/GSEA_CONvCG_0.5.GseaPreranked.1578923228842/gsea_report_for_na_pos_1578923228842.xls",sep="\t")
gsea_con_v_cg0.5_pos <- gsea_con_v_cg0.5_pos[,-c(2,3,12)]
colnames(gsea_con_v_cg0.5_pos) <- unlist(gsea_con_v_cg0.5_pos[1,])
gsea_con_v_cg0.5_pos <- gsea_con_v_cg0.5_pos[-1,]

gsea_con_v_cg0.5_neg <- read.table("GSEA/GSEA_CONvCG_0.5/GSEA_CONvCG_0.5.GseaPreranked.1578923228842/gsea_report_for_na_neg_1578923228842.xls",sep="\t")
gsea_con_v_cg0.5_neg <- gsea_con_v_cg0.5_neg[,-c(2,3,12)]
colnames(gsea_con_v_cg0.5_neg) <- unlist(gsea_con_v_cg0.5_neg[1,])
gsea_con_v_cg0.5_neg <- gsea_con_v_cg0.5_neg[-1,]

gsea_con_v_cg0.5 <- rbind(gsea_con_v_cg0.5_pos,gsea_con_v_cg0.5_neg)

# Control vs CG 1
gsea_con_v_cg1.0_pos <- read.table("GSEA/GSEA_CONvCG_1.0/GSEA_CONvCG_1.0.GseaPreranked.1578923228839/gsea_report_for_na_pos_1578923228839.xls",sep="\t")
gsea_con_v_cg1.0_pos <- gsea_con_v_cg1.0_pos[,-c(2,3,12)]
colnames(gsea_con_v_cg1.0_pos) <- unlist(gsea_con_v_cg1.0_pos[1,])
gsea_con_v_cg1.0_pos <- gsea_con_v_cg1.0_pos[-1,]

gsea_con_v_cg1.0_neg <- read.table("GSEA/GSEA_CONvCG_1.0/GSEA_CONvCG_1.0.GseaPreranked.1578923228839/gsea_report_for_na_neg_1578923228839.xls",sep="\t")
gsea_con_v_cg1.0_neg <- gsea_con_v_cg1.0_neg[,-c(2,3,12)]
colnames(gsea_con_v_cg1.0_neg) <- unlist(gsea_con_v_cg1.0_neg[1,])
gsea_con_v_cg1.0_neg <- gsea_con_v_cg1.0_neg[-1,]

gsea_con_v_cg1.0 <- rbind(gsea_con_v_cg1.0_pos,gsea_con_v_cg1.0_neg)

# Collect NES scores for all comparisons
df1 <- gsea_con_v_mg1.5[,c("NAME","NES")]
colnames(df1) <- c("path_id","NES_mg1.5")
df2 <- gsea_con_v_mg3.0[,c("NAME","NES")]
colnames(df2) <- c("path_id","NES_mg3.0")
df3 <- gsea_con_v_mg6.0[,c("NAME","NES")]
colnames(df3) <- c("path_id","NES_mg6.0")
df4 <- gsea_con_v_cg0.5[,c("NAME","NES")]
colnames(df4) <- c("path_id","NES_cg0.5")
df5 <- gsea_con_v_cg1.0[,c("NAME","NES")]
colnames(df5) <- c("path_id","NES_cg1.0")

# Merge NES scores from different comparisons into one data frame
gsea_nes <- merge(merge(merge(merge(df1,df2,by="path_id"),df3,by="path_id"),df4,by="path_id"),df5,by="path_id")
gsea_nes$path_id <- tolower(gsea_nes$path_id)
rownames(gsea_nes) <- gsea_nes$path_id
gsea_nes <- gsea_nes[,-1]

# Collect FDR scores for all comparisons
df1 <- gsea_con_v_mg1.5[,c("NAME","FDR q-val")]
colnames(df1) <- c("path_id","FDR_mg1.5")
df2 <- gsea_con_v_mg3.0[,c("NAME","FDR q-val")]
colnames(df2) <- c("path_id","FDR_mg3.0")
df3 <- gsea_con_v_mg6.0[,c("NAME","FDR q-val")]
colnames(df3) <- c("path_id","FDR_mg6.0")
df4 <- gsea_con_v_cg0.5[,c("NAME","FDR q-val")]
colnames(df4) <- c("path_id","FDR_cg0.5")
df5 <- gsea_con_v_cg1.0[,c("NAME","FDR q-val")]
colnames(df5) <- c("path_id","FDR_cg1.0")

# Merge FDR scores from different comparisons into one data frame
gsea_fdr <- merge(merge(merge(merge(df1,df2,by="path_id"),df3,by="path_id"),df4,by="path_id"),df5,by="path_id")
gsea_fdr$path_id <- tolower(gsea_fdr$path_id)
rownames(gsea_fdr) <- gsea_fdr$path_id
gsea_fdr <- gsea_fdr[,-1]

# Resulting data frame of FDRs and NESs consist of factors, convert factors to numeric
indx <- sapply(gsea_fdr, is.factor)
gsea_fdr[indx] <- lapply(gsea_fdr[indx], function(x) as.numeric(as.character(x)))
indx <- sapply(gsea_nes, is.factor)
gsea_nes[indx] <- lapply(gsea_nes[indx], function(x) as.numeric(as.character(x)))
# Reverse the sign of NES to perform the comparison to control
gsea_nes <- -gsea_nes

# Transform FDR for plotting heatmap
fdr_pos_nes <- 1-gsea_fdr
fdr_neg_nes <- -1+gsea_fdr

# Create a matrix for plotting a heatmap 
fdr_heatmap <- matrix(,nrow=nrow(gsea_fdr),ncol=ncol(gsea_fdr))
for (i in seq(nrow(gsea_fdr))){
  for (j in seq(ncol(gsea_fdr))){
    # If NES >0, then the element of the matrix equals 1-FDR
    if (gsea_nes[i,j]>0){
      fdr_heatmap[i,j] <- fdr_pos_nes[i,j]
    } else { # If NES <0, then the element of the matrix equals -1+FDR
      fdr_heatmap[i,j] <- fdr_neg_nes[i,j]
    }
  }
}

fdr_heatmap <- data.frame(fdr_heatmap)

rownames(fdr_heatmap) <- rownames(gsea_fdr)
colnames(fdr_heatmap) <- colnames(gsea_fdr)

# Get pathway names
fdr_heatmap$path_id <- rownames(fdr_heatmap)
fdr_heatmap <- merge(path_id_name_human,fdr_heatmap,by="path_id")
rownames(fdr_heatmap) <- fdr_heatmap$path_name
fdr_heatmap <- fdr_heatmap[,-c(1,2)]
fdr_heatmap <- fdr_heatmap[,c(4,5,1,2,3)]

# Remove pathways that have FDR >=0.05 for all comparisons
fdr_heatmap <- fdr_heatmap[which(rowSums(abs(fdr_heatmap)>0.95)>=1),]

hc <- hclust(d = dist(fdr_heatmap, method = "euclidean"))

# Plot dendogram based on clustering results
dendo <- as.dendrogram(hc)

# Plot dendro
dendro.plot <- ggdendrogram(data = dendo, rotate = TRUE)
print(dendro.plot)

fdr_heatmap <- fdr_heatmap[hc$order,]

# Compose a matrix for plotting heatmap
com_heat <- expand.grid(path_id = rownames(fdr_heatmap),cond = colnames(fdr_heatmap))
# Add variable: transformed FDR
com_heat$fdr <- as.vector(as.matrix(fdr_heatmap))

p <- (ggplot(data=com_heat, aes(x = cond, y = path_id)) 
  + geom_tile(aes(fill = fdr))
  # indicate paramters of the colorbar of the legend
  + scale_fill_gradient2(low = "blue", high = "red", mid = "white",midpoint = 0, space = "Lab",name="fdr")
  # put the legend on the top 
  + theme(legend.position="top",
          legend.text=element_text(family="Helvetica",color="black", size=8),
          axis.title.x=element_blank(),axis.title.y=element_blank(),
          axis.text.x = element_text(family="Helvetica",size=8,angle=45,hjust = 1),
          axis.text.y = element_text(family="Helvetica",color="black", size=8))
  # put gene labels on the right and reverse the order
  + scale_y_discrete(position="right",labels=rownames(fdr_heatmap))
  + scale_x_discrete(labels=c("Collagen (0.5 mg/ml)","Collagen (1.0 mg/ml)","Matrigel (1.5 mg/ml)","Matrigel (3.0 mg/ml)","Matrigel (6.0 mg/ml)"))
  + coord_equal())

print(p)










































