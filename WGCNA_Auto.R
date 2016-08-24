#display current working direcotry
getwd()

#If necessary, give path to where data is stored
workingDir = "/Users/rebeccakrasnoff/Documents/Current/Willsey/POGZ_Eirene/Data/pogz_P2/"
setwd(workingDir);

#Load WGCNA package (if not already loaded)
library(WGCNA)

#following is important - do not omit
options(stringsAsFactors = FALSE)

#Allow multi-threading within WGCNA. Helps speed up certain calculations.
#At present this call is necessary for the code to work.
#Any error here may be ignored, but you may want to update WGCNA if you see one.
#Caution: skip this line if you run Rstudio or other third-party R envirnments
#enableWGCNAThreads()

#Load the data saved in the first script
lnames = load(file = "pogz_P2-01-dataInput.RData")
#the variable lnames contains the names of loaded variables
lnames

#####################################################################################################
#Automatic construction of the gene network and identification of modules
#Choosing the soft-thresholding power: analysis of network topology
#choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from=12, to=20, by=2))
#call the network topology analysis function
sft = pickSoftThreshold(dataExpr, powerVector=powers, verbose=5)
#plot the results
sizeGrWindow(6,5)
par(mfrow = c(1,2))
cex1 = 0.9
#Scale-free topolgy fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
                               xlab="Soft Thresdhold (power)", ylab = "Scale Free Topology Model Fit, signed R^2",
                               type = "n", main = paste("Scale indeendence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
     labels=powers, cex=cex1, col="red")
#this line corresponds to using an R^2 cut off of h
abline(h=0.30, col="red") ###look into this number some more....

##################################################################################################
#One step network construction and module detection
net = blockwiseModules(datExpr, power = 6, 
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25, 
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE, saveTOMFileBase = "pogz_P2_TOM",
                       verbose = 3)
###Word of caution - mess with other values for this function; many have been left
#as default, which may not be optimal for other data

#To see how many modules were identified and what module sizes are
table(net$colors)

#display dendrogram with color assignent
#open a graphics window
sizeGrWindow(10,5)
dev.new()
#convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
#Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], 
                    "Module colors", dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

##################################################################################
#save module assignment and module eigengene information necessary for futher analysis
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs
geneTree = net$dendrograms[[1]]
workingDir = "/Users/rebeccakrasnoff/Documents/Current/Willsey/POGZ_Eirene/Data/pogz_P2/"
setwd(workingDir);
save(MEs, moduleLabels, moduleColors, geneTree, 
     file = "pogz_P2_networkConstruction_auto.RData")


###########################################################################################



