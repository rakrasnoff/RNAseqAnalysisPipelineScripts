options(stringsAsFactors = FALSE)

wrkdir <- "/Users/rebeccakrasnoff/Documents/Current/Willsey/Hypoxia/Data/times_merged/genes.fpkm_tracking"

myGenes <- read.delim(wrkdir, stringsAsFactors = F)
head(myGenes)

myIndex <- apply(myGenes[,grep("status", colnames(myGenes))], 1, function(x) {
    status <- grepl("LOWDATA", x)
    !any(status)
    })
passingGenes <- myGenes[myIndex, ]
nrow(myGenes)
nrow(passingGenes)
head(passingGenes)

#get rid of unusable genes: mir, RP##, AC, commas
useGenes <- passingGenes[grep("*mir*", passingGenes$gene_short_name, invert=TRUE),]
nrow(useGenes)
useGenes <- useGenes[grep("^RP.*-", useGenes$gene_short_name, invert=TRUE),]
nrow(useGenes)
useGenes <- useGenes[grep(",", useGenes$gene_short_name, invert=TRUE),]
nrow(useGenes)
useGenes <- useGenes[grep("^AC", useGenes$gene_short_name, invert=TRUE),]
nrow(useGenes)
#are still duplicate gene names - can get rid of if necessary
#useGenes <- unique(useGenes)

passingGenesTrimmed <- useGenes[, grep("FPKM", colnames(passingGenes))]
geneIDS <- useGenes$gene_short_name
length(passingGenesTrimmed)
nrow(passingGenesTrimmed)

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

