#GSA Lightening

install.packages("devtools")
library("devtools")
install_github("billyhw/GSALightning")

library(GSALightning)
source("https://bioconductor.org/biocLite.R")
biocLite("ELMER")

data(df.glite)
data(sampleInfo)
data(targetGenes)
#Remove genes with 0 sample variance
expression <- expression[apply(expression, 1, sd) !=0,]

#Run GSALight

GSALightResults <- GSALight(eset=df.glite, fac = factor)



