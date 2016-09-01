#RNA seq analysis pipeline
#steps:
#1: load cufflinks data, 
#2: load cuffnorm data
#3: clean data
#4: isolate differentially expressed genes
#5: convert gene ids
#6: GO analysis
#7: GSA Light
#8: WGCNA analysis
#to add eventually: load Hseq input

#code pieces

#choose project
projectName <- readline("What project are you working on? (Current options: Eirene, or POGZ)")  
chooseProject <- function(projectName) {
   if (projectName == "Eirene" | projectName == "POGZ" | projectName == "pogz") {
    setwd("/Users/rebeccakrasnoff/Documents/Current/Willsey/POGZ_Eirene/Data")
  }
}


folderName <- readline("What condition do you want to run for? 
                       (Current options: p2, e16)") 
chooseFile <- function(folderName){
  if (foldername == "p2") {
    setwd("pogz_P2")
    outputPath <- "/Users/rebeccakrasnoff/Documents/Current/Willsey/POGZ_Eirene/Results/pogz_P2"
  }
  if (folderName == "e16") {
    setwd("pogz_e16")
    outputPath <- "/Users/rebeccakrasnoff/Documents/Current/Willsey/POGZ_Eirene/Results/pogz_e16"
  }
}

#reset values (if running again)
loadCufflinks = 0
loadCuffnorm = 0
cleanSeqData = 0
findDiff = 0
convertIds = 0
runGO = 0
runGSA = 0


chooseOptions <- function(){
  optionNames <- readline("What analysis do you want to run? 
(1) Load cufflinks data
(2) Load cuffnorm data
(3) Clean data (or, define background) 
(4) Find differentially expressed genes
(5) Convert gene ids
(6) GO analysis
(7) GSA analysis") 
  if ("1" %in% optionNames) {
   loadCufflinks = 1
  }
  if ("2" %in% optionNames) {
    loadCuffnorm = 1
  }
  if ("3" %in% optionNames) {
    cleanSeqData = 1
  }
  if ("4" %in% optionNames) {
    findDiff = 1
  }
  if ("5" %in% optionNames) {
    convertIds = 1
  }
  if ("6" %in% optionNames) {
    runGo = 1
  }
  if ("7" %in% optionNames) {
    runGSA = 1
  }
}

setwd("/Users/rebeccakrasnoff/Documents/Current/Willsey/RNAseqAnalysisPipeline/RNAseqAnalysisPipelineScripts/")

if (loadCufflinks) {
  source("R/loadCufflinks.R")
} 



fun <- function(){
  x <- readline("What is the value of x?")  
  y <- readline("What is the value of y?")
  t <- readline("What are the T values?")
  v <- readline("What are the V values?")
  
  #x <- as.numeric(unlist(strsplit(x, ",")))
  #y <- as.numeric(unlist(strsplit(y, ",")))
  #t <- as.numeric(unlist(strsplit(t, ",")))
  #v <- as.numeric(unlist(strsplit(v, ",")))
  
  out1 <- x + y
  out2 <- t + v
  
  return(list(out1, out2))
  
}




