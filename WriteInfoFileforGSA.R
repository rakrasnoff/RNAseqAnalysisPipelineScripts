#Create file with sample info for cuffnorm data to be used in GSALite

pogzP2info <- data.frame("sample" = c("WT1", "WT2", "WT3", "HET1", "HET2", "HET3")) #this row can be anything, can just be 1:n where n is number of samples
pogzP2info$type <- as.factor(c("Control", "Control", "Control", "Experiment", "Experiment", "Experiment"))
#set wd to set file location
setwd("/Users/rebeccakrasnoff/Documents/Current/Willsey/POGZ_Eirene/Data/pogz_P2/")
write.table(pogzP2info, file = "P2_info_file.txt", sep = '/t', row.names=FALSE, col.names=TRUE, quote = FALSE)
#for more samples, use:
control.vec <- rep("Control", n/2)
experimentl.vec <- rep("Experiment", n/2)


#creat gene set
genes_mouse <- read.table("/Users/rebeccakrasnoff/Documents/Current/Willsey/RNAseqAnalysisPipeline/Mouse_to_human_gene_conversionKey.csv", sep=",", header=TRUE)
asd_genes <- read.table("~/Documents/Current/Willsey/HEK293expression/Data/TADA_top65ASDgenes_June2016.csv", sep=",", header=TRUE)
top_reg <- read.table("~/Documents/Current/Willsey/HEK293expression/Data/TADA_top179ASDgenes_July2016_annotatedByFunction.csv", sep=",", header=TRUE)

head(genes_mouse)
head(asd_genes)
head(top_reg)

reg_mouse <- reg_mouse$Mouse
asd_mouse <- asd_mouse$Mouse


head(gene_list)
gene_list_t <- t(gene_list)
head(gene_list_t)

##
gene_sets = list(top_65=asd_mouse, top_reg=reg_mouse)






