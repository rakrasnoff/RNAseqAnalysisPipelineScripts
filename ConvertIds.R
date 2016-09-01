source("http://bioconductor.org/biocLite.R")
biocLite(org.Mm.eg.db)
#biocLite('org.Hs.eg.db')
library(org.Hs.eg.db)
library(org.Mm.eg.db)

#Convert Gene IDs
wrkdir = ""
setwd(wrkdir)
load("cuffLinksLoaded.RData")

#Convert differentially expressed gene ids
#Give List of Genes as Vector

#show some genes
head(gene.diff$gene)

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
head(select(org.Mm.eg.db, as.vector(gene.diff$gene), "ENTREZID", "SYMBOL"),50)

#create vector of converted ids
convertedDiff<- select(org.Mm.eg.db, as.vector(gene.diff$gene), "ENTREZID", "SYMBOL")

#append vector to dataframe if applicable
#list.diff$entrez_id <- convertedids$ENTREZID[match(gene.exp$gene, convertedids$SYMBOL)]

#write file
write.table(convertedDiff, file="convertedDiff.txt", quote=FALSE, sep = "\t", row.names = FALSE)


#Convert background gene ids

head(gene.exp$gene)

##### Look at list of possible id types
keytypes(org.Mm.eg.db)
#current.type <- "SYMBOL"
#convert.type <- "ENTREZID"
#convert ids
head(select(org.Mm.eg.db, as.vector(gene.exp$gene), "ENTREZID", "SYMBOL"),10)

#create vector of converted ids
convertedExp<- select(org.Mm.eg.db, as.vector(gene.exp$gene), "ENTREZID", "SYMBOL")

#append vector to dataframe if applicable
#list.diff$entrez_id <- convertedids$ENTREZID[match(gene.exp$gene, convertedids$SYMBOL)]

#write file
write.table(convertedExp, file="convertedExp.txt", quote=FALSE, sep = "\t", row.names = FALSE)

###For only up regulated genes

head(gene.up$gene)
#current.type <- "SYMBOL"
#convert.type <- "ENTREZID"
head(select(org.Mm.eg.db, as.vector(gene.up$gene), "ENTREZID", "SYMBOL"),10)

#create vector of converted ids
convertedUp<- select(org.Mm.eg.db, as.vector(gene.up$gene), "ENTREZID", "SYMBOL")

#write file
write.table(convertedUp, file="convertedUp.txt", quote=FALSE, sep = "\t", row.names = FALSE)

##for only down regulated genes
head(gene.down$gene)
#current.type <- "SYMBOL"
#convert.type <- "ENTREZID"
head(select(org.Mm.eg.db, as.vector(gene.down$gene), "ENTREZID", "SYMBOL"),5)

#create vector of converted ids
convertedDown<- select(org.Mm.eg.db, as.vector(gene.down$gene), "ENTREZID", "SYMBOL")

#write file
write.table(convertedDown, file="convertedDown.txt", quote=FALSE, sep = "\t", row.names = FALSE)


##################################
##Save data
save(convertedDiff, convertedExp, convertedUp, convertedDown, file = "convertedIds.RData")
