source("https://bioconductor.org/biocLite.R")
biocLite("ReactomePA")
biocLite("clusterProfiler")

library(ReactomePA)
#data(geneList)
#head(geneList)
#length(geneList)
#nrow(geneList)
library(org.Hs.eg.db)
library(org.Mm.eg.db)

#load data from id-converting-script
getwd()
wrkdir <- "/Users/rebeccakrasnoff/Documents/Current/Willsey/POGZ_Eirene/Data/pogz_P2"
setwd(wrkdir)

load("convertedIds.RDATA")

universe <- convertedExp$ENTREZID
gene <- convertedDiff$ENTREZID

x <- enrichPathway(as.vector(gene), organism = "mouse", pvalueCutoff = 0.05,
                   pAdjustMethod = "fdr", qvalueCutoff = 0.2, as.vector(universe), minGSSize = 10,
                   maxGSSize = 500, readable = FALSE)
head(summary(x))
results <- x

#visualize results
barplot(x, showCategory = 8)
dotplot(x, showCategory = 15)
enrichMap(x, layout=igraph::layout.kamada.kawai, vertex.label.cex = 1)
#cnetplot(x, categorySize="pvalue", foldChange=geneList)

#compare clusters
require(clusterProfiler)
res <- compareCluster(gene, fun="enrichPathway")
plot(res)

