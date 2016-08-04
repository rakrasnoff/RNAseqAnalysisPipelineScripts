########## Load required packages ##########
library(gplots)
library(ggplot2)
library(gridExtra)
library(WGCNA)
library(impute)
library(GO.db)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(annotate)
library(compare)
library(cummeRbund)

####### Set working directory
setwd("/Users/rebeccakrasnoff/Documents/Current/Willsey/POGZProject/data/cuffdiff_out") ##write path to file here

##background
gene.exp <- gene.exp.full[
  (gene.exp.full$value_1 >= 1 | gene.exp.full$value_2 >= 1)
  & gene.exp.full$status == "OK"
  & gene.exp.full$gene != "-", ]
nrow(gene.exp)
nrow(gene.exp.full)
##### Define foreground ("gene list") and background genes
#Parameters for FG genes: difference in expression was >+2x, q value was less than .05,
#Note: if not done above, also filter out genes without names, with expression values less than one both sample 1 and 2, and with status !OK
gene.FG <- gene.exp[
  (gene.exp$log2.fold_change. >= 1 | gene.exp$log2.fold_change. <= -1) 
  & gene.exp$q_value < .05, ]
head(gene.FG)
#BG genes: using whole list of genes with names, expression values >1, status OK
gene.BG <- gene.exp

#### Gene lists are ready for functional annotation
### FUNCTIONAL ANNOTATION: UP-REGULATED AND DOWN-REGULATED GENES 

dyn.load('/Library/Java/JavaVirtualMachines/jdk1.8.0_91.jdk/Contents/Home/jre/lib/server/libjvm.dylib')
require(rJava)
.jinit()
.jcall("java/lang/System", "S", "getProperty", "java.runtime.version")
library("RDAVIDWebService")

###Connect to webservice
# Create a DAVIDWebService object connected to David, using your registration email.
# To register, go to: http://david.abcc.ncifcrf.gov/content.jsp?file=WS.html.
david <- DAVIDWebService(email="willseylab@ucsf.edu", url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")


####RUN UPREGULATED AND DOWN REGULATED GENES THROUGH DAVID SIMULTANEOUSLY
# Define foreground and background gene lists.
# The foreground list should be contained within the background list.
# note: for possible inputs for idType, use getIdTypes(david)
FG <- addList(david, gene.FG$entrez_id, idType="ENTREZ_GENE_ID", listName="foreground", listType="Gene")
BG <- addList(david, gene.BG$entrez_id, idType="ENTREZ_GENE_ID", listName="background", listType="Background")

# Inspect FG and BG to see the proportion of genes recognized by David, and those that are unmapped.
FG 
BG

# Inspect "david" object to see the gene lists selected as foreground and background.
david

# Specifiy annotation categories.
setAnnotationCategories(david, c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL"))

# Get functional annotation chart as R object.
FuncAnnotChart <- getFunctionalAnnotationChart(david)
View(FuncAnnotChart)
# Print functional annotation chart to file.
getFunctionalAnnotationChartFile(david, "FuncAnnotChart.tsv")

# Get functional annotation clustering (limited to 3000 genes).
FuncAnnotClust <- getClusterReport(david)
head(summary(FuncAnnotClust))
dev.new()
plot2D(FuncAnnotClust, 1)

# Print functional annotation clustering to file (limited to 3000 genes).
getClusterReportFile(david, "FuncAnnotClust.tsv")

###############
davidGODag<-DAVIDGODag(members(FuncAnnotClust)[[1]], pvalueCutoff=0.05, "CC")
plotGOTermGraph(g=goDag(davidGODag),r=davidGODag, max.nchar=40, node.shape="ellipse")
dev.new()
davidGODag<-DAVIDGODag(members(FuncAnnotClust)[[2]], pvalueCutoff=0.1, "CC")
plotGOTermGraph(g=goDag(davidGODag),r=davidGODag, max.nchar=40, node.shape="ellipse")


####RUN ONLY UP REGULATED GENES THROUGH DAVID 
gene.u <- gene.FG[
  (gene.FG$log2.fold_change. >= 1) , ]
nrow(gene.u)
FGu <- addList(david, gene.u$entrez_id, idType="ENTREZ_GENE_ID", listName="upregulated", listType="Gene")
BG <- addList(david, gene.BG$entrez_id, idType="ENTREZ_GENE_ID", listName="background", listType="Background") #this line is only here because it wouldn't let me set the same backgroun as above
david

# Specifiy annotation categories.
setAnnotationCategories(david, c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL"))

# Get functional annotation chart as R object.
FuncAnnotChart.up <- getFunctionalAnnotationChart(david)
View(FuncAnnotChart.up)
nterms <- .1*(nrow(FuncAnnotChart.up))
FuncAnnotGraph.up <- data.frame(FuncAnnotChart.up$Term[1:nterms])
FuncAnnotGraph.up[,2] <- data.frame(FuncAnnotChart.up$PValue[1:nterms])
#PLOTTTTTTT

# Print functional annotation chart to file.
getFunctionalAnnotationChartFile(david, "FuncAnnotChart_up.tsv")

# Get functional annotation clustering (limited to 3000 genes).
FuncAnnotClust <- getClusterReport(david)
head(summary(FuncAnnotClust))
dev.new()
plot2D(FuncAnnotClust, 1)
davidGODag<-DAVIDGODag(members(FuncAnnotClust)[[1]], pvalueCutoff=0.05, "CC")
plotGOTermGraph(g=goDag(davidGODag),r=davidGODag, max.nchar=40, node.shape="ellipse")


####RUN ONLY DOWN REGULATED GENES THROUGH DAVID 
gene.d <- gene.FG[
  (gene.FG$log2.fold_change. <= (-1)) , ]
nrow(gene.d)
FGd <- addList(david, gene.d$entrez_id, idType="ENTREZ_GENE_ID", listName="downregulated", listType="Gene")
BG <- addList(david, gene.BG$entrez_id, idType="ENTREZ_GENE_ID", listName="backgroundu", listType="Background")
david


#set categories
setAnnotationCategories(david, c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL"))

# Get functional annotation chart as R object.
FuncAnnotChart.d <- getFunctionalAnnotationChart(david)
View(FuncAnnotChart.d)
nterms <- .1*(nrow(FuncAnnotChart.d))
FuncAnnotGraph.d <- data.frame(FuncAnnotChart.d$Term[1:nterms])
FuncAnnotGraph.d[,2] <- data.frame(FuncAnnotChart.d$PValue[1:nterms])
##pick up here

# Print functional annotation chart to file.
getFunctionalAnnotationChartFile(david, "FuncAnnotChart_d.tsv")

# Get functional annotation clustering (limited to 3000 genes).
FuncAnnotClust <- getClusterReport(david)
head(summary(FuncAnnotClust))
dev.new()
plot2D(FuncAnnotClust, 1)
davidGODag<-DAVIDGODag(members(FuncAnnotClust)[[1]], pvalueCutoff=0.05, "CC")
plotGOTermGraph(g=goDag(davidGODag),r=davidGODag, max.nchar=40, node.shape="ellipse")