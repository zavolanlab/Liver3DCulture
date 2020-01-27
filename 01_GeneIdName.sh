# Download human GTF
wget ftp://ftp.ensembl.org/pub/release-96/gtf/homo_sapiens/Homo_sapiens.GRCh38.96.gtf.gz
gunzip Homo_sapiens.GRCh38.96.gtf.gz

# Get genes associated with chromosomes
awk '($1 == "1") || ($1 == "2") || ($1 == "3") || ($1 == "4") || ($1 == "5") || ($1 == "6") || ($1 == "7") || ($1 == "8") || ($1 == "9") || ($1 == "10") || ($1 == "11") || ($1 == "12") || ($1 == "13") || ($1 == "14") || ($1 == "15") || ($1 == "16") || ($1 == "17") || ($1 == "18") || ($1 == "19") || ($1 == "20")|| ($1 == "21") || ($1 == "22") || ($1 == "23") || ($1 == "MT") || ($1 == "X")  || ($1 == "Y") { print $0}' Homo_sapiens.GRCh38.96.gtf > tmp2
awk '($3 == "gene") { print $10"\t"$14}' tmp2 > tmp3

# Assign Ensembl ids to gene symbols
sed 's/\"//g' tmp3 > tmp4
sed 's/;//g' tmp4 > GeneIdTOGeneName.txt

# Remove temporary variables
rm tmp2
rm tmp3
rm tmp4
