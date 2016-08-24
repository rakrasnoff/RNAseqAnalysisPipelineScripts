# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored.
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = "/Users/rebeccakrasnoff/Downloads/FemaleLiver-Data/"
setwd(workingDir);
# Load the WGCNA package
library(WGCNA);
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

#Read in my dataset - using dataset from load cuffnorm
#femData = read.csv("LiverFemale3600.csv") - rather than loading in, using genes.norm for P2 right now
# Take a quick look at what is in the data set:
dim(genes.norm)
names(genes.norm)
fix(genes.norm)

#The expression data set contains 6 samples. 
#Note that each row corresponds to a gene and column to a sample or auxiliary information. 
#We now transpose the expression data for further analysis.
datExpr0 = as.data.frame(t(genes.norm));
names(datExpr0) = rownames(genes.norm)
rownames(datExpr0) = names(genes.norm)
dim(datExpr0)

# We first check for genes and samples with too many missing values:
gsg = goodSamplesGenes(datExpr0, verbose=3)
gsg$allOK

# if gsg$allOK is not true:

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

#Now, cluster samples to see if there are any obvious outliers
sampleTree = hclust(dist(datExpr0), method = "average")
#Plot the sample tree: open a graphic output window of size 12 by 9 inches
#User should change dimsnions if window is too large or too small
#I changed it to 6x5
sizeGrWindow(6,5)
#pdf(file="Plots/sampleClustering.pdf", width = 6, height = 5)
par(cex = .6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

#choose threshold to remove outlier - #HET1 was outlier
#plot line to show cut
abline(h=200000, col="red")
#determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 200000, minSize = 0)
table(clust)
#clust 1 contains the samples we want to keep 
###Because I have so few samples, this next bit was not helpful, instead, just manually remove outlier
#keepSamples = (clust == 1)
#datExpr = datExpr0[keepSamples, ]
#nGenes = ncol(datExpr)
#nSamples = nrow(datExpr)

#manually remove outlier
rownames(datExpr0)
datExpr = datExpr0[-4, ]
rownames(datExpr)
#the variable datExpr now contains the expression data ready for network analysis

#Load clinical trait data
#read the trait data and match the samples for which they were measured to the expression samples
#

traitData = read.delim("/Users/rebeccakrasnoff/Documents/Current/Willsey/POGZ_Eirene/Data/pogz_P2/P2_info_file.txt")
rownames(traitData) = traitData$sample
dim(traitData)
names(traitData)

#remove columns with info we don't need - already done for mine, so allTraits just equals traitData
allTraits = traitData
rownames(allTraits) = traitData$sample
dim(allTraits)
names(allTraits)

#create a data frame analogous to expression data that will hold the clinical traits
sampleNames = rownames(datExpr)
traitRows = match(sampleNames, allTraits$sample)
datTraits = allTraits[traitRows,] #this gets rid of any samples removed from the cluster in my trait dat
rownames(datTraits) = allTraits[traitRows, 1]

collectGarbage()

#Now have expression data stored in datExpr, and clinical traits in datTraits. Now,
#visualize how clinical traits relate to sample dendrogram

#recluster samples
sampleTree2 = hclust(dist(datExpr), method="average")
#convert traits to color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(as.numeric(datTraits$type), signed=FALSE)
#plot the sample dendrogram and colors underneath
plotDendroAndColors(sampleTree2, traitColors, groupLabels=names(datTraits), 
                    main = "Sample Dendrogram with Trait Heatmap")

#this dendrograph essentially tells me that my HET are experiment samples, and my control are control samples

#Now, save relevant expression and trait data for use in the next steps of this tutorial
workingDir = "/Users/rebeccakrasnoff/Documents/Current/Willsey/POGZ_Eirene/Data/pogz_P2/"
setwd(workingDir);
save(datExpr, datTraits, file = "pogz_P2-01-dataInput.RData")


