#load cuff diff for shiny app

gene.exp.full <- read.table(input$filename, header=TRUE)
min_fpkm = 1 #specify minimum value for fpkm for samples 1 and 2; for no minimum, enter 0
status_ok = "OK" #filter out by status value? If not, change to ... (figure out best way to do this)
require_gene_name = "-" #will only load genes with official gene symbol

gene.exp <- gene.exp.full[
  (gene.exp.full$value_1 >= input$min_fpkm | gene.exp.full$value_2 >= input$min_fpkm)
  & gene.exp.full$gene != "-",
  , ]
if (status_ok == TRUE) {
  gene.exp[gene.exp$status == "OK")
}

nrow(gene.exp)
nrow(gene.exp.full)