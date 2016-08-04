##trying different formats to see if I can get Gsea r to work

#Load in this data just to create dummy list
HETfpkm1 <- read.table("/Users/rebeccakrasnoff/Documents/Current/Willsey/POGZ_EireneData/Data/cufflinks_out/pogzP2-het1_ACTTGA_L001_001/genes.fpkm_tracking", header=TRUE)
HETfpkm2 <- read.table("/Users/rebeccakrasnoff/Documents/Current/Willsey/POGZ_EireneData/Data/cufflinks_out/pogzP2-het2_GATCAG_L001_001/genes.fpkm_tracking", header=TRUE)
HETfpkm3 <- read.table("/Users/rebeccakrasnoff/Documents/Current/Willsey/POGZ_EireneData/Data/cufflinks_out/pogzP2-het3_GGCTAC_L001_001/genes.fpkm_tracking", header=TRUE)
WTfpkm1 <- read.table("/Users/rebeccakrasnoff/Documents/Current/Willsey/POGZ_EireneData/Data/cufflinks_out/pogzP2-wt1_ATCACG_L001_001/genes.fpkm_tracking", header=TRUE)
WTfpkm2 <- read.table("/Users/rebeccakrasnoff/Documents/Current/Willsey/POGZ_EireneData/Data/cufflinks_out/pogzP2-wt2_TTAGGC_L001_001/genes.fpkm_tracking", header=TRUE)
WTfpkm3 <- read.table("/Users/rebeccakrasnoff/Documents/Current/Willsey/POGZ_EireneData/Data/cufflinks_out/pogzP2-wt3_TAGCTT_L001_001/genes.fpkm_tracking", header=TRUE)




NAME = HETfpkm1$gene_id #give column with gene name #these are just the cuff ids for now
Description <- HETfpkm1$locus #give column with description of gene or gene symbol - locus for now

#Create table
df.gsea.fake <- data.frame("WT1"= WTfpkm1$FPKM[1:25302], "WT2"= WTfpkm2$FPKM[1:25302], 
                           "WT3" = WTfpkm3$FPKM[1:25302], "HET1" = HETfpkm1$FPKM, 
                           "HET2" = HETfpkm2$FPKM[1:25302], "HET3" = HETfpkm3$FPKM[1:25302], row.names = HETfpkm1$locus)
head(df.gsea.fake)

write.table(df.gsea.fake, file="gsalitefake.csv", quote = FALSE, row.names=TRUE)





######create fake gct withougt 0s

head(df.gsea.fake.nonzero)
nrow(df.gsea.fake.nonzero)

### try without 0s
n <- nrow(df.gsea.fake.nonzero)
cl <- length(df.gsea.fake.nonzero)
n.samples <- 6 #number of sample groups
setwd("/Users/rebeccakrasnoff/Documents/Current/Willsey/Scripts/GSEAPR/POGZ_dataset")

#write gct
write.table("#1.2", file="GSEA_dummy.gct", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(t((c(n, n.samples))), file="GSEA_dummy.gct", sep="\t", row.names=FALSE, col.names=FALSE, append=TRUE, quote=FALSE)
write.table(df.gsea.fake.nonzero, file="GSEA_dummy.gct", sep="\t", row.names=FALSE, col.names=TRUE, append=TRUE, quote=FALSE)

#####create cls file
write.table(t((c(n.samples, 2, 1))), file="GSEA_dummy.cls", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
sample1.name <- "WT"
sample2.name <- "HET"
write.table(t((c("#", sample1.name, sample2.name))), file="GSEA_dummy.cls", sep="\t", row.names=FALSE, col.names=FALSE, append=TRUE, quote=FALSE)
sample1.list <- rep("WT", (n.samples/2))
sample2.list <- rep("HET", (n.samples/2))
write.table(t((c(sample1.list, sample2.list))), file="GSEA_dummy.cls", sep="\t", row.names=FALSE, col.names=FALSE, append=TRUE, quote=FALSE)


#look at some of their data files to compare




## trying different formats to see if I can get GSA Lite to work


