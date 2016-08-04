#Load Data From Cufflinks, Filter data

#location of your cuff diff file
setwd("/Users/rebeccakrasnoff/Documents/Current/Willsey/POGZProject/data/cuffdiff_out") ##write path to file here

#Load data in
gene.exp.full <- read.table("gene_exp.diff", header=TRUE)
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

#write new data

write.table(gene.exp, file="gene.exp.csv", header=TRUE, quote=FALSE, sep = ",")
