setwd("~/Dropbox/WillseyLab/collaborations/willsey_helen/xenopus_ASD/gene_expression/")

# source("http://bioconductor.org/biocLite.R")
# biocLite(c("AnnotationDbi", "impute", "GO.db", "preprocessCore"))
# biocLite("org.Hs.eg.db")
# install.packages("WGCNA")

library(ggplot2)
library(WGCNA)
library(reshape2)

options(stringsAsFactors = FALSE)
allowWGCNAThreads()
################################################################################
#LOAD RAW DATA
################################################################################

rawData <- read.delim(file = "owensPaper/clutchA_polya_absolute_TPE_gene_isoform.txt", 
                      header = T)
tmpData <- subset(rawData, X!="Isoform")
tmpData <- tmpData[ , -1]

myTime <- unlist(tmpData[1 ,  -1])
myGeneIDs <- tmpData[-1 , 1]

tmpData <- tmpData[-1 , -1]
rownames(tmpData) <- myGeneIDs

myGeneSymbols <- sapply(myGeneIDs, function(x) strsplit(x, "|", fixed = T)[[1]][2]) 

tmpData <- tmpData[!is.na(myGeneSymbols) & myGeneSymbols != "unnamed", ]
myGeneSymbols <- myGeneSymbols[!is.na(myGeneSymbols) & myGeneSymbols != "unnamed"]

################################################################################
#QC RAW DATA
################################################################################

nas <- apply(tmpData, 1, is.na)
any(nas) #FALSE

range(tmpData)
summary(tmpData)

#perhaps think about function to remove outlier gene expression values (i.e. outside 95% CI or 3 SDs)?
cleanData <- tmpData
rm(list = c("tmpData"))

################################################################################
#WGCNA
################################################################################

load("/Users/rebeccakrasnoff/Documents/Current/Willsey/Hypoxia/Data/pasca-hypoxia-dataInput.RData")

# save
# save(datExpr, myTime, myGeneIDs, file = "ClutchA_polyA_TPE_Absolute_dataInput.RData")


##############
# module construction
# 
#load("ClutchA_polyA_TPE_Absolute_dataInput.RData")

#using datExpr from load dat wgcna script


# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#calculate adjacency matrix using chosen soft threshold power
softPower = 3;
adjacency = adjacency(datExpr, power = softPower);

# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);

# want large modules, so set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

# Convert numeric labels into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

#merge similar modules based on correlation between eigengenes
# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

#choose height cut of 0.25, corresponding to correlation of 0.75 for merging
MEDissThres = 0.05
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

#check merging graphically
sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#dev.off()

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, 
     file = "pasca-hypoxia-networkConstruction.RData")

load("pasca-hypoxia-networkConstruction.RData")


toPlot <- MEs
# myTime <- myTime[-which(names(myTime)=="Sample_76")]
toPlot$time <- myTime
# toPlot$row <- rep(1:3, each=30)[1:nrow(toPlot)]

toPlot <- melt(toPlot, id.vars = c("time"), value.name="expression", variable.name = "module")
toPlot$color <- gsub(pattern = "ME", replacement = "", x = toPlot$module)
toPlot$row <- rep(1:3, each=89)

p <- ggplot(toPlot, aes(time, expression)) + geom_point(col=toPlot$color) + geom_line()
p + facet_wrap(row ~ module) + theme_bw()

# ggsave(p, height=10, width=30)

#group genes by module

myModuleGenes <- sapply(colnames(datExpr0), function(x) strsplit(x, "|", fixed = T)[[1]][2])
myModuleGenes <- myModuleGenes[!is.na(myModuleGenes)]

table(moduleColors)

myModuleGenes <- data.frame(gene <- myModuleGenes, module <- moduleColors)
names(myModuleGenes) <- c("gene", "module")
myModuleGenes$XenBase <- NA
#need to load myGeneKey
myModuleGenes$XenBase[myModuleGenes$gene %in% myGeneKey$GeneSymbol] <- myGeneKey$XB_Gene[match(myModuleGenes$gene, myGeneKey$GeneSymbol, nomatch = 0)]

#check modules for ASD gene enrichment

#load TADA list
tada <- read.delim(file = "~/Dropbox/WillseyLab/geneFiles/TADA_SNV_CNV_combined_Feb7.txt", header = T)
asdGene <- subset(tada, qvalue.combined < 0.1, select = RefSeqName)
asdGeneWithOrtholog <- subset(myModuleGenes, gene %in% tolower(asdGene$RefSeqName))
names(asdGeneWithOrtholog) <- c("gene", "module")

#calculate enrichment by module using hypergeometric

myModules <- data.frame(module <- unique(moduleColors))
colnames(myModules) <- c("module")
myModules$size <- table(moduleColors)[match(myModules$module, names(table(moduleColors)))]
myModules$asdGenes <- 0

myModules$asdGenes[myModules$module %in% names(table(asdGeneWithOrtholog$module))] <- 
  table(asdGeneWithOrtholog$module)[match(myModules$module, names(table(asdGeneWithOrtholog$module)), nomatch = 0 )] 

myModules$enrichmentPval <- apply(myModules[,2:3], 1, function(x) phyper(x[2] - 1, 65, 14405 - 65, x[1], lower.tail = F) )

asdModules <- subset(myModules, asdGenes>0)

write.table(subset(myModuleGenes, module=="bisque4", select="gene"), file = "bisque4_DAVID.txt", 
            row.names = F, col.names = F, quote = F)
write.table(subset(myModuleGenes, module=="lightyellow", select="XenBase"), file = "lightyellow_DAVID.txt", 
            row.names = F, col.names = F, quote = F)
write.table(subset(myModuleGenes, module=="skyblue3", select="gene"), file = "skyblue3_DAVID.txt", 
            row.names = F, col.names = F, quote = F)
write.table(subset(myModuleGenes, module=="lightcyan", select="gene"), file = "lightcyan_DAVID.txt", 
            row.names = F, col.names = F, quote = F)

toPlot <- toPlot[toPlot$color %in% asdModules$module,]

p <- ggplot(toPlot, aes(time, value)) + geom_point(col=toPlot$color) + geom_line() #+ annotate(geom = "text", x=20, y=0.25, label = toPlot)
p + facet_grid(. ~ variable)

asdGene$RefSeqName[which(tolower(asdGene$RefSeqName) %in% myModuleGenes$gene[myModuleGenes$module=="lightyellow"])]
asdGene$RefSeqName[which(tolower(asdGene$RefSeqName) %in% myModuleGenes$gene[myModuleGenes$module=="bisque4"])]
asdGene$RefSeqName[which(tolower(asdGene$RefSeqName) %in% myModuleGenes$gene[myModuleGenes$module=="lightcyan"])]
asdGene$RefSeqName[which(tolower(asdGene$RefSeqName) %in% myModuleGenes$gene[myModuleGenes$module=="brown4"])]
asdGene$RefSeqName[which(tolower(asdGene$RefSeqName) %in% myModuleGenes$gene[myModuleGenes$module=="skyblue3"])]