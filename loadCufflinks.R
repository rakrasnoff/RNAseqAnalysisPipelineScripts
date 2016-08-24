#choose file based on filename
loadCufflinks <- function(projectName, fileName) {
  if (projectName == "Eirene" | projectName == "POGZ" | projectName == "pogz") {
    setwd("/Users/rebeccakrasnoff/Documents/Current/Willsey/POGZ_Eirene/Data")
  }
  if (foldername == "p2") {
    setwd("pogz_P2")
  }
  if (folderName == "e16") {
    setwd("pogz_e16")
  }
  cufflinksData <- read.table("gene_exp.diff", header=TRUE)
}



