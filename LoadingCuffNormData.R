#Script to easily load cuffnorm output into environment. Set up to be customizable and easy to use

#path to cuffnorm genes.fpkm_tracking file
#give sample name (P2, or e16)
sample = "P2"

if (sample == "P2") {
cuff.norm <- read.table("/Users/rebeccakrasnoff/Documents/Current/Willsey/POGZ_Eirene/Data/pogz_P2/cuffnorm_out/genes.fpkm_tracking", header=TRUE)
}
if (sample == "e16") {
cuff.norm <- read.table("/Users/rebeccakrasnoff/Documents/Current/Willsey/POGZ_Eirene/Data/pogz_e16/genes.fpkm_tracking", header=TRUE)
}

#look at data
head(cuff.norm)
#get column names for data
colnames(cuff.norm)


#Redefine column names so they are easier to access, trunkate columns in the dataframe with necessary information
genes.norm <- data.frame("gene.name" = cuff.norm$gene_short_name)
genes.norm$gene.id <- cuff.norm$gene_id
genes.norm$WT1 <- cuff.norm[,10]
genes.norm$WT1.status <- cuff.norm[,13]
genes.norm$WT2 <- cuff.norm[,14]
genes.norm$WT2.status <- cuff.norm[,17]
genes.norm$WT3 <- cuff.norm[,18]
genes.norm$WT3.status <- cuff.norm[,21]
genes.norm$HET3 <- cuff.norm[,22]
genes.norm$HET3.status <- cuff.norm[,25]
genes.norm$HET2 <- cuff.norm[,26]
genes.norm$HET2.status <- cuff.norm[,29]
genes.norm$HET1 <- cuff.norm[,30]
genes.norm$HET1.status <- cuff.norm[,33]
nrow(genes.norm)
head(genes.norm)


#filter where status is not OK
#filter where name doesn't exist
genes.norm <- genes.norm[genes.norm$gene.name != "-" & genes.norm$WT1.status == "OK" & genes.norm$WT2.status == "OK"
                         & genes.norm$WT3.status == "OK" & genes.norm$HET1.status == "OK" & genes.norm$HET2.status == "OK"
                         & genes.norm$HET3.status == "OK",]
nrow(genes.norm)

#look for ASD names in duplicates
nrow(genes.norm)
nrow(genes.norm[duplicated(genes.norm$gene.name),])
nrow(genes.norm[unique(genes.norm$gene.name),])
duplicates <- data.frame(genes.norm$gene.name[duplicated(genes.norm$gene.name)])
nrow(unique(duplicates))
asd_genes <- read.table("~/Documents/Current/Willsey/HEK293expression/Data/TADA_top65ASDgenes_June2016.csv", sep=",", header=TRUE)
head(asd_genes)
asd_genes$RefSeqName %in% duplicates

#reformat dataframe with gene names as row names
head(genes.norm)
nrow(genes.norm)
gene.name.v <- as.character(genes.norm$gene.name)
length(gene.name.v)
gene.name.u <- make.unique(gene.name.v, sep = ".")
length(gene.name.u)
genes.norm <- data.frame("WT1" = genes.norm$WT1, "WT2" = genes.norm$WT2, "WT3" = genes.norm$WT3, "HET1" = genes.norm$HET1, "HET2" = genes.norm$HET2, "HET3" = genes.norm$HET3, row.names = gene.name.u)
head(genes.norm) 

#write to file

if (sample == "P2") {
  write.table(genes.norm, file="/Users/rebeccakrasnoff/Documents/Current/Willsey/POGZ_Eirene/Data/pogz_P2/pogzP2GSA.txt", quote=FALSE, sep = "\t", row.names = TRUE)
  
}

if (sample == "e16") {
  write.table(genes.norm, file="/Users/rebeccakrasnoff/Documents/Current/Willsey/POGZ_Eirene/Data/pogz_e16/pogze16GSA.txt", quote=FALSE, sep = "\t", row.names = TRUE)
}

}


