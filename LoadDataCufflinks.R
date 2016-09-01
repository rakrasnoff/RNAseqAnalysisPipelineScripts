#Load Data From Cufflinks, Filter data

#location of your cuff diff file

#Load data in
wrkdir <- "/Users/rebeccakrasnoff/Documents/Current/Willsey/POGZ_Eirene/Data/pogz_P2/"
setwd(wrkdir)
gene.exp.full <- read.delim("gene_exp.diff")
head(gene.exp.full)

#Sort Data 
#DAVID script will sort more, this is just to clean up the data

min_fpkm = 1 #specify minimum value for fpkm for samples 1 and 2; for no minimum, enter 0
status_value = "OK" #filter out by status value? If not, change to ... (figure out best way to do this)
require_gene_name = "-" #will only load genes with official gene symbol

gene.exp <- gene.exp.full[
  (gene.exp.full$value_1 >= min_fpkm | gene.exp.full$value_2 >= min_fpkm)
  & gene.exp.full$status == status_value
  & gene.exp.full$gene != require_gene_name, ]

nrow(gene.exp)
nrow(gene.exp.full)

#gene.exp is now a cleaned dataset
print("gene.exp is now a clean dataset")

#identify differentially expressed genes
gene.diff <- gene.exp[
  (gene.exp$log2.fold_change. >= 1 | gene.exp$log2.fold_change. <= -1) 
  & gene.exp$q_value < .05, ]
list.diff <- gene.diff$gene

nrow(gene.diff)

#write new data

write.table(gene.exp, file="geneExp.txt", quote=FALSE, sep = "\t", row.names = FALSE)
write.table(gene.diff, file="geneDiff.txt", quote=FALSE, sep="\t", row.names = FALSE)
write.table(list.diff, file = "listDiff.txt", quote=FALSE, sep="\t", row.names = FALSE)

## which genes are upregulated?
gene.up <- gene.diff[
  (gene.diff$log2.fold_change. >= 1) , ]
nrow(gene.up)

## which genes are downregulated?

gene.down <- gene.diff[
  (gene.diff$log2.fold_change. <= -1) , ]
nrow(gene.down)

#write to tables
write.table(gene.up, file="geneUp.txt", quote=FALSE, sep = "\t", row.names = FALSE)
write.table(gene.down, file="geneDown.txt", quote=FALSE, sep = "\t", row.names = FALSE)
save(gene.exp, gene.diff, gene.up, gene.down, file = "cuffLinksLoaded.RData")

