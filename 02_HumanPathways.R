#!/usr/bin/env Rscript

### Created: Jan 27, 2020
### Author: Anastasiya Boersch
### Company: Zavolan Group, Biozentrum, University of Basel

# Load required libraries
source("https://bioconductor.org/biocLite.R")
biocLite("KEGGREST")
library(KEGGREST) # The version KEGGREST_1.16.1 was used
library(biomaRt)

path.list <- keggLink("hsa","pathway")
path.list.m <- as.matrix(path.list)
path.id <- unique(rownames(path.list.m))
path.gene <- data.frame(rownames(path.list.m),path.list.m)
colnames(path.gene) <- c("path_id", "entrez_id")
ent.id <- substring(path.gene$entrez_id,5)
path.gene$entrez_id <- as.numeric(ent.id)

# Convert Entrez Ids to Ensembl gene Ids
mart = useMart('ensembl')
listDatasets(mart)
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
x <- listAttributes(human)
mapping <- getBM(attributes = c("entrezgene_id", "ensembl_gene_id"), mart = human, uniqueRows = TRUE)
colnames(mapping) <- c("entrez_id","ensembl_id")
map_ens_uniq <- mapping[!duplicated(mapping$ensembl_id), ]

# Merge pathways, entrez ids and pathway ids
path.gene <- merge(path.gene,map_ens_uniq,by="entrez_id")

# Add gene names
ens2name <- read.table("GeneIdTOGeneName.txt",header=F, sep="\t")
colnames(ens2name) <- c("ensembl_id","gene_symbol")

path.gene <- merge(path.gene,ens2name,by="ensembl_id")

# Add pathway names
pathid <- as.character(unique(path.gene$path_id))
tmp1 <- rep(NA,length(pathid))
tmp2 <- rep(NA,length(pathid))
for (i in seq(length(pathid))) {
  tmp2[i] <- pathid[i]
  res <- try(keggGet(pathid[i])[[1]]$NAME)
  if(!inherits(res, "try-error"))
  { tmp1[i] <- keggGet(pathid[i])[[1]]$NAME
  tmp1[i] <- substr(tmp1[i],1,nchar(tmp1[i])-23)
  }
}
path_id_name <- data.frame(path_id=tmp2,path_name=tmp1)

# Write a table containing pathway IDs and pathway names
write.table(path_id_name,"HumanPathIdName.txt",sep="\t",quote = T)

# Write a table containing pathways and genes annotated to these pathways
path.gene <- merge(path.gene,path_id_name,by="path_id")
write.table(path.gene,"HumanGenesPathways.txt",sep="\t",quote = T)
