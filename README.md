This repository contains all scripts used for RNA-Seq data analysis presented in the manuscript by Ghosh et al. "The transcriptional landscapes of a hepatoma cell line grown on scaffolds of extracellular matrix proteins".

To get started, clone the repository and change to the directory you specified as the working directory:
```bash
git clone https://github.com/zavolanlab/Liver3DCulture.git path/to/workdir
cd path/to/workdir
```

Download and unzip pre-processed RNA-Seq data "GSE142206_counts.txt.gz" (a text file containing raw counts per gene across conditions) and "GSE142206_tpms.txt.gz" (a text file containing normalized expression per gene in TPM units across conditions) published in GEO in the working directory (see https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE142206). The analysis is based on processing these data files.  

Execute the following scripts in the indicated order:
1. Create a table matching ensembl gene ids and gene symbols by running:
```bash
sh 01_GeneIdGeneName.sh
```
The resulting table is saved as "GeneIdTOGeneName.txt".
2. Create a list of human KEGG pathways and associated genes. Run R script "02_HumanPathways.R". Resulting tables are saved as "HumanGenesPathways.txt" and "HumanPathIdName.txt".     
3. Perform differential expression (DE) analysis. Run R script "03_DE_analysis.R". A directory with results called "DE_analysis" will be created within your working directory.  
4. Perform gene set enrichment analysis. Run R script "04_GSEA.R". A directory with GSEA results called "GSEA" will be created within your working directory. Note that before running the script you should install GSEA software for the command line, see https://software.broadinstitute.org/gsea/downloads.jsp  
5. Perform principal component analysis. Run R script "05_PCA.R".