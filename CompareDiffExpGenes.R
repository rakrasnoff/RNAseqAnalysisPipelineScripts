#compare differentially expressed genes across samples

######Load dataset 1
gene.exp.full <- read.table("/Users/rebeccakrasnoff/Documents/Current/Willsey/POGZ_Eirene/Data/Pogz_P2/gene_exp.diff", header=TRUE)
head(gene.exp.full)

#Sort Data 

min_fpkm = 1 #specify minimum value for fpkm for samples 1 and 2; for no minimum, enter 0
status_value = "OK" #filter out by status value? If not, change to ... (figure out best way to do this)
require_gene_name = "-" #will only load genes with official gene symbol

gene.exp1 <- gene.exp.full[
  (gene.exp.full$value_1 >= min_fpkm | gene.exp.full$value_2 >= min_fpkm)
  & gene.exp.full$status == status_value
  & gene.exp.full$gene != require_gene_name, ]

nrow(gene.exp1)
nrow(gene.exp.full)

#gene.exp is now a cleaned dataset

#identify differentially expressed genes
gene.diff1 <- gene.exp1[
  (gene.exp1$log2.fold_change. >= 1 | gene.exp1$log2.fold_change. <= -1) 
  & gene.exp1$q_value < .05, ]
list.diff1 <- gene.diff1$gene

nrow(gene.diff1)


#######Load dataset 2
gene.exp.full <- read.table("/Users/rebeccakrasnoff/Documents/Current/Willsey/POGZ_Eirene/Data/pogz_e16/gene_exp.diff", header=TRUE)
head(gene.exp.full)

#Sort Data 
#DAVID script will sort more, this is just to clean up the data

min_fpkm = 1 #specify minimum value for fpkm for samples 1 and 2; for no minimum, enter 0
status_value = "OK" #filter out by status value? If not, change to ... (figure out best way to do this)
require_gene_name = "-" #will only load genes with official gene symbol

gene.exp2 <- gene.exp.full[
  (gene.exp.full$value_1 >= min_fpkm | gene.exp.full$value_2 >= min_fpkm)
  & gene.exp.full$status == status_value
  & gene.exp.full$gene != require_gene_name, ]

nrow(gene.exp2)
nrow(gene.exp.full)

#gene.exp is now a cleaned dataset

#identify differentially expressed genes
gene.diff2 <- gene.exp2[
  (gene.exp2$log2.fold_change. >= 1 | gene.exp2$log2.fold_change. <= -1) 
  & gene.exp2$q_value < .05, ]
list.diff2 <- gene.diff2$gene

nrow(gene.diff2)

#compare lists
overlap <- intersect(list.diff1, list.diff2)

write.table(overlap, file="/Users/rebeccakrasnoff/Documents/Current/Willsey/POGZ_Eirene/Data/pogz_e16/GenesDiffExpInP2and16.csv", row.names = FALSE, quote=FALSE)


## create lists of up regulated genes

gene.u1 <- gene.diff1[
  (gene.diff1$log2.fold_change. >= 1) , ]
nrow(gene.u1)

gene.u2 <- gene.diff2[
  (gene.diff2$log2.fold_change. >= 1) , ]
nrow(gene.u2)

#list of down regulated genes 

gene.d1 <- gene.diff1[
  (gene.diff1$log2.fold_change. <= (-1)) , ]
nrow(gene.d1)

gene.d2 <- gene.diff2[
  (gene.diff2$log2.fold_change. <= (-1)) , ]
nrow(gene.d2)

#compare up regulated genes
overlap.up <- intersect(gene.u1$gene, gene.u2$gene)

#compare down regulated genes
overlap.down <- intersect(gene.d1$gene, gene.d2$gene)

intersect(gene.u1$gene, gene.d2$gene)
e16Up2D <- intersect(gene.u2$gene, gene.d1$gene)


#save data
save(e16Up2D, overlap.up, overlap.down, overlap, file = "/Users/rebeccakrasnoff/Documents/Current/Willsey/POGZ_Eirene/Data/sampleComparison/OverlappingData.RData")
#write tables
setwd("/Users/rebeccakrasnoff/Documents/Current/Willsey/POGZ_Eirene/Data/pogz_e16/")
write.table(overlap.up, file="/Users/rebeccakrasnoff/Documents/Current/Willsey/POGZ_Eirene/Data/pogz_e16/GenesUpExpInP2and16.csv", row.names = FALSE, quote=FALSE)
write.table(overlap.down, file="/Users/rebeccakrasnoff/Documents/Current/Willsey/POGZ_Eirene/Data/pogz_e16/GenesDownExpInP2and16.csv", row.names = FALSE, quote=FALSE)


#convert Ids
library(org.Hs.eg.db)
library(org.Mm.eg.db)

#Convert Gene IDs
wrkdir = ""
setwd(wrkdir)
load("OverlappingData.RData")

#choose list to convert:
toConvert <- overlap
#Set configuration
##What type are your current ids?
current.type <- "SYMBOL"
##What kind of ids do you want to convert to?
convert.type <- "ENTREZID"
#convert ids
#create vector of converted ids
conv_overlap<- select(org.Mm.eg.db, as.vector(toConvert), "ENTREZID", "SYMBOL")
#for others
toConvert <- overlap.up #none here
conv_overlapUp <- select(org.Mm.eg.db, as.vector(toConvert), "ENTREZID", "SYMBOL")
#
toConvert <- overlap.down
conv_overlapDown <- select(org.Mm.eg.db, as.vector(toConvert), "ENTREZID", "SYMBOL")
#
toConvert <- e16Up2D
conv_e16Up2D <- select(org.Mm.eg.db, as.vector(toConvert), "ENTREZID", "SYMBOL")

save(conv_overlap, conv_overlapDown, conv_e16Up2D, file = "/Users/rebeccakrasnoff/Documents/Current/Willsey/POGZ_Eirene/Data/sampleComparison/Conv_OverlappingData.RData")
