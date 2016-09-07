options(stringsAsFactors = FALSE)
        
#load data
wrkdir <- "/Users/rebeccakrasnoff/Dropbox/Hypoxia_Data/Data/times_merged/genes.fpkm_tracking"
myGenes <- read.delim(wrkdir, stringsAsFactors = F)

#look at data
head(myGenes)
#right now, genes are rows, and samples are columns (sample_FPKM  sample_conf_lo  sample_conf_hi  sample_status)

############DATA CLEANING##############
#create new data frame with only genes where status is "OK" for all samples
myIndex <- apply(myGenes[,grep("status", colnames(myGenes))], 1, function(x) {
  status <- grep("OK", x, invert=TRUE)
  !any(status)
})
okGenes <- myGenes[myIndex, ]
nrow(myGenes)
nrow(okGenes)
head(okGenes)
#now, data only has genes where status is okay for all samples

##For sencing tables to J
#write.table(myGenes$gene_short_name, file = "/Users/rebeccakrasnoff/Desktop/allGeneNamesHypoxia.txt", sep = "\t", quote = F)
#write.table(myGenes, file = "/Users/rebeccakrasnoff/Desktop/allGenesHypoxia.txt", , sep = "\t", quote = F)
#duplicates <- subset(myGenes, duplicated(gene_short_name))
#write.table(duplicates$gene_short_name, file = "/Users/rebeccakrasnoff/Desktop/dupNamesGenesHypoxia.txt", sep = "\t", quote = F)
#write.table(duplicates, file = "/Users/rebeccakrasnoff/Desktop/dupsGenesHypoxia.txt", sep = "\t", quote = F)

#Get rid of rows that are not genes: 
gGenes <- okGenes
badGeneNames <- data.frame("*mir*", "-", "^RP.*-", ",", "^AC", "RNA", "sno", "U8", "U7", "U3", "U6", "U1", "_SRP", "SNORD")
getGenes <- function(x) {
  gGenes <<- data.frame(gGenes[grep(x, gGenes$gene_short_name, invert=TRUE),])
  print(nrow(gGenes))
}
apply(badGeneNames, 2, getGenes)

#useGenes <- subset(useGenes, !duplicated(gene_short_name)) #this would get rid of duplicates

###create new dataset with only FPKMS
#data frame containing only gene Id (first column) and fpkm values for all samples
geneIdFPKM <- data.frame("geneId" = gGenes$gene_short_name, gGenes[, grep("FPKM", colnames(gGenes))])

#####
#get rid of genes with <.2 coefficient of variation and/or no FPKM of at least 2 in at least 1 sample
passingGenes <- geneIdFPKM[apply(geneIdFPKM[,-1],1,var) >= .2 & apply(geneIdFPKM[,-1],1,max)>=2,]
genesWGCNA <- passingGenes[, grep("FPKM", colnames(passingGenes))]
geneIds <- passingGenes$geneId
#format for GSA
passingGenes$geneId[duplicated(passingGenes$geneId)]
genes.gsa <- subset(passingGenes, !duplicated(passingGenes$geneId))
genes.gsa <- data.frame(genes.gsa[,-1], row.names = genes.gsa$geneId)
head(genes.gsa)

#load sample info to match trait data
metadata <- read.delim("/Users/rebeccakrasnoff/Documents/Current/Willsey/Hypoxia/Data/times_merged/orderedSampleDetails.txt", stringsAsFactors = F, header = F)
metadata <- unique(metadata)

#create dataframe of trait data
sampleIds <- colnames(passingGenesTrimmed)
sampleNums <- do.call(rbind, strsplit(sampleIds, "_"))
sampleNumIds <- cbind(sampleIds, sampleNums)
sampleNumIds <- data.frame(sampleNumIds[,c(1,3)])
names(sampleNumIds) <- c("sampleName", "sampleNum")
traitDat0 <- merge(sampleNumIds, metadata, by.x = "sampleNum", by.y = "V1")
rownames(traitDat0) <- traitDat0$sampleName
traitDat <- traitDat0[,c(2,3,5)]


datExpr0 <- t(passingGenesTrimmed)

#save
save(datExpr0, traitDat, geneIDS, file = "pasca-hypoxia-cuffnorm-to-WGCNA.RData")



#####

library(data.table)