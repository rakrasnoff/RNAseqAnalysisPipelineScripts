#WGCNA analysis of top-co-expressed genes
#############################################################################################################################################

#get list of genes best correlated to each index gene
getTopCor<-function(indexGene, expressionData, topNumber){
  #calculate pairwise correlations between index gene and all other genes
  index=(which(colnames(datExpr) %in% indexGene))
  tmp=cor(expressionData[,index], datExpr, method="pearson")
  #tmp[lower.tri(tmp, diag=TRUE)]<-NA
  #tmp[abs(tmp)<=0.7]<-NA
  #order by correlation, choose top20 that are also abs(correlation)>0.7
  topCorr=order(abs(tmp), decreasing=T)[2:(topNumber+1)]
  geneCor=tmp[topCorr]
  names(geneCor)=colnames(tmp)[topCorr]
  geneCor=geneCor[abs(geneCor)>=0.7]
  names(geneCor)=
    return(geneCor)
}

#using list of network genes, get all connections >0.7
getNetworkConnections<-function(networkGenes, expressionData){
  #calculate pairwise correlations between index gene and all other genes
  index=(which(colnames(datExpr) %in% networkGenes))
  tmp=cor(expressionData[,index], use="pairwise.complete.obs", method="pearson")
  tmp[lower.tri(tmp, diag=TRUE)]<-NA
  tmp[abs(tmp)<=0.7]<-NA
  index2=which(!is.na(tmp), arr.ind=T)
  pairs=cbind(colnames(tmp)[index2[,2]],tmp[index2], rownames(tmp)[index2[,1]])
  return(pairs)
}

#genes dex in the same direction versus control across 5 or more conditions
sharedDex=c("1810008I18Rik","2010003K11Rik","Arntl","C330021F23Rik","Ccl2","Cdkn1a","Cidec","Cyp26b1","Cyp2b10",
            "Foxq1","Gadd45b","Gm16551","Hsd3b5","Hspa1b","Klhdc7a","Mycn","Nlrp12","Pdzk1ip1","Relb","Srebf1","Tuba8","Ube2c", "Chrna4")

myTopCorr=lapply(sharedDex, function(x) getTopCor(x, datExpr, 10))
myBestGenes=unique(names(unlist(myTopCorr)))
genePairs=getNetworkConnections(myBestGenes, datExpr)

write.table(genePairs, "testGenepairs.txt", row.names=F, col.names=F, append=F, quote=F, sep="\t")
#############################################################################################################################################



# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
#plot results
setEPS()
postscript(file=paste(myStamp, "choosePower.eps", sep="_"), width=10, height=5, paper="special", colormodel="srgb", horizontal=FALSE)
par( ps=12, font.lab=2, font.main=2, omi=c(0.5,0,0,0.5), mfrow = c(1,2),
     mgp=c(2.75, 0.5, 0), las=0, cex=1, cex.lab=1, cex.main=1, cex.axis=1)
cex1=0.9
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
dev.off()

#choose 6 as softthresholding power (SFT.R.sq=0.941), then calculate adjacencies
softPower = 6
#unsigned
# adjacency = adjacency(datExpr, power = softPower, type="unsigned")
#signed
adjacency = adjacency(datExpr, power = softPower, type="signed")
# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM

#use hierarchical clustering to create dendrogram of gene-gene relationships
# Call the hierarchical clustering function
geneTree = flashClust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
# plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
# labels = FALSE, hang = 0.04)
#identify modules by dynamic branch cutting
# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize ...