source("http://bioconductor.org/biocLite.R")
biocLite(org.Mm.eg.db)
#biocLite('org.Hs.eg.db')
library(org.Hs.eg.db)
library(org.Mm.eg.db)

#Convert Gene IDs

#Give List of Genes as Vector
gene_list <- gene.exp$gene

#show some genes
head(gene_list)

##### Look at list of possible id types
keytypes(org.Mm.eg.db)

#Set configuration
##What type are your current ids?
current.type <- "SYMBOL"
##What kind of ids do you want to convert to?
convert.type <- "ENTREZID"

#make sure you have selected the correct ids for current and convert
head(keys(org.Mm.eg.db, keytype=current.type))
head(keys(org.Mm.eg.db, keytype=convert.type))

#convert ids
head(select(org.Mm.eg.db, as.vector(gene_list), "ENTREZID", "SYMBOL"),50)

#create vector of converted ids
convertedids<- select(org.Mm.eg.db, as.vector(gene_list), "ENTREZID", "SYMBOL")

#append vector to dataframe if applicable
gene.exp$entrez_id <- convertedids$ENTREZID[match(gene.exp$gene, convertedids$SYMBOL)]

#write file
write.table(gene.exp, file="/Users/rebeccakrasnoff/Documents/Current/Willsey/POGZ_Eirene/Data/pogz_e16/geneExpConv.csv", quote=FALSE, sep = ",", row.names = FALSE)

