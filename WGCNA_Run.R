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








