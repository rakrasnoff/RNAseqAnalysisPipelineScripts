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
#Read in the female liver data set
femData = read.csv("LiverFemale3600.csv")
# Take a quick look at what is in the data set:
dim(femData)
names(femData)
fix(femData)

#The expression data set contains 135 samples. 
#Note that each row corresponds to a gene and column to a sample or auxiliary information. 
#We now remove the auxiliary data and transpose the expression data for further analysis.
datExpr0 = as.data.frame(t(femData[, -c(1:8)]));
names(datExpr0) = femData$substanceBXH
rownames(datExpr0) = names(femData)[-c(1:8)]
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

#choose threshold to remove outlier
#plot line to show cut
abline(h=15, col="red")
#determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 15, minSize = 10)
table(clust)
#clust 1 contains the samples we want to keep
keepSamples = (clust == 1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

#the variable datExpr now contains the expression data ready for network analysis

#Load clinical trait data
#read the trait data and match the samples for which they were measured to the expression samples

traitData = read.csv("ClinicalTraits.csv")
dim(traitData)
names(traitData)

#remove columns with info we don't need
allTraits = traitData[, -c(31, 16)]
allTraits = allTraits[, c(2, 11:36) ]
dim(allTraits)
names(allTraits)

#create a data frame analogous to expression data that will hld the clinical traits
femaleSamples = rownames(datExpr)
traitRows = match(femaleSamples, allTraits$Mice)
datTraits = allTraits[traitRows, -1]
rownames(datTraits) = allTraits[traitRows, 1]

collectGarbage()

#Now have expression data stored in datExpr, and clinical traits in datTraits. Now,
#visualize how clinical traits relate to sample dendrogram

#recluster samples
sampleTree2 = hclust(dist(datExpr), method="average")
#convert traits to color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed=FALSE)
#plot the sample dendrogram and colors underneath
plotDendroAndColors(sampleTree2, traitColors, groupLabels=names(datTraits), 
                    main = "Sample Dendrogram with Trait Heatmap")

#Now, save relevant expression and trait data for use in the next steps of this tutorial
save(datExpr, datTraits, file = "FemaleLiver-01-dataInput.RData")


