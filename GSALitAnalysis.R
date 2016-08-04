#GSA Lightning Analysis

#tutorial
library(devtools) 
#install_github("billyhw/GSALightning")
library(GSALightning)
library(ELMER)
library(ELMER.data)
data(expression)
data(sampleInfo)
data(targetGenes)

expression <- expression[apply(expression,1,sd) != 0,]

GSALightResults <- GSALight(eset = expression, fac = factor(sampleInfo$TN), gs = targetGenes, 
                            nperm = 1000, minsize = 10, rmGSGenes = 'gene')
head(GSALightResults)

#with my data
gene.exp <- read.table("/Users/rebeccakrasnoff/Documents/Current/Willsey/POGZ_EireneData/Data/cuffdiff_out/gene_exp.diff", header=TRUE)
head(fakeData)
gene.exp <- gene.exp[gene.exp$status == "OK"
  & gene.exp$gene != "-", ]

gene.exp <- gene.exp[gene.exp$gene == unique(gene.exp$gene),]

list1 <- gene.exp$gene == unique(gene.exp$gene)

uniquegenes <- unique(gene.exp$gene)

gene.gsa <- data.frame("WT" = gene.exp$value_1, "HET" = gene.exp$value_2, row.names = gene.exp$gene)
nrow(gene.exp)
nrow(gene.exp.full)

fakeData <- fakeData[apply(fakeData,1,sd) != 0,]
GSALightResults <- GSALight(eset = expression, fac = factor(sampleInfo$TN), gs = targetGenes, 
                            nperm = 1000, minsize = 10, rmGSGenes = 'gene')
