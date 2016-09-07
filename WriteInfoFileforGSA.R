#Create file with sample info for cuffnorm data to be used in GSALite

pogzP2info <- data.frame("sample" = c("WT1", "WT2", "WT3", "HET1", "HET2", "HET3")) #this row can be anything, can just be 1:n where n is number of samples
pogzP2info$type <- as.factor(c("Control", "Control", "Control", "Experiment", "Experiment", "Experiment"))
#set wd to set file location
setwd("/Users/rebeccakrasnoff/Documents/Current/Willsey/POGZ_Eirene/Data/pogz_P2/")
write.table(pogzP2info, file = "P2_info_file.txt", sep = '\t', row.names=FALSE, col.names=TRUE, quote = FALSE)
#for more samples, use:
control.vec <- rep("Control", n/2)
experimentl.vec <- rep("Experiment", n/2)


#gs
genes_mouse <- read.table("/Users/rebeccakrasnoff/Documents/Current/Willsey/RNAseqAnalysisPipeline/Mouse_to_human_gene_conversionKey.csv", sep=",", header=TRUE)
asd_genes <- read.table("~/Documents/Current/Willsey/HEK293expression/Data/TADA_top65ASDgenes_June2016.csv", sep=",", header=TRUE)
top_reg <- read.table("~/Documents/Current/Willsey/HEK293expression/Data/TADA_top179ASDgenes_July2016_annotatedByFunction.csv", sep=",", header=TRUE)

head(genes_mouse)
head(asd_genes)
head(top_reg)


#for asd genes
asd_mouse <- merge(genes_mouse, asd_genes, by.x="Human", by.y="RefSeqName")
reg_mouse <- merge(genes_mouse, top_reg, by.x="Human", by.y="ID")

write.table(asd_mouse, file="/Users/rebeccakrasnoff/Documents/Current/Willsey/RNAseqAnalysisPipeline/RNAseqAnalysisPipelineScripts/asd_mouse.txt", quote=FALSE, sep="\t", row.names=F, col.names=FALSE)
write.table(reg_mouse, file="/Users/rebeccakrasnoff/Documents/Current/Willsey/RNAseqAnalysisPipeline/RNAseqAnalysisPipelineScripts/reg_mouse.txt", quote=FALSE, sep="\t", row.names=F, col.names=FALSE)

head(reg_mouse)
head(asd_mouse)
reg_mouse <- reg_mouse$Mouse
asd_mouse <- asd_mouse$Mouse


##
gene_sets = list(top_65=asd_mouse, top_reg=reg_mouse)




#write.table(gene_sets, file="/Users/rebeccakrasnoff/Documents/Current/Willsey/RNAseqAnalysisPipeline/GSA_Lightning.textClippingGSAgenesets.csv", quote=FALSE, sep=",", row.names=TRUE, col.names=FALSE)


