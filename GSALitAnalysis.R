#GSA Lightning Analysis


library(devtools) 
#install_github("billyhw/GSALightning")
library(GSALightning)
library(data.table)
library(ELMER)
library(ELMER.data)
##for tutorial
#data(expression)
#data(sampleInfo)
#data(targetGenes)


#load necessary files
gene_sets <- gene_sets = list(top_65=asd_mouse, top_reg=reg_mouse)
info_file <- read.delim("/Users/rebeccakrasnoff/Documents/Current/Willsey/POGZ_Eirene/Data/pogz_P2/P2_info_file.txt")
#genes.norm <- read.table("/Users/rebeccakrasnoff/Documents/Current/Willsey/POGZ_Eirene/Data/pogz_P2/pogzP2GSA.txt")
#genes.norm <- read.table("/Users/rebeccakrasnoff/Documents/Current/Willsey/POGZ_Eirene/Data/pogz_e16/pogze16GSA.txt")
genes.gsa <- genes.norm[apply(genes.norm,1,sd) != 0,]

GSALigthtResults <- GSALight(eset = genes.gsa, fac = factor(info_file$type), gs = gene_sets, nperm = 10000,
                            minsize = 1, rmGSGenes = 'gene', verbose=TRUE)


head(GSALightResults)

setwd("/Users/rebeccakrasnoff/Documents/Willsey/POGZ_Eirene/Data/pogz_P2")

write.table(GSALightResults, file='P2GSALightResults.txt', sep = "\t", quote = FALSE)

