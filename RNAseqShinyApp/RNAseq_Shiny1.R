library(shiny)

ui <- fluidPage(
  titlePanel("RNAseq Analysis"),
  #create side bar
  sidebarLayout(
    sidebarPanel(
  fileInput(inputId = "filename", label = "Load file here:", multiple = FALSE, 
            accept = NULL, width = NULL),
  numericInput(inputId = "n_sample", label = "How many samples does your data have?", value = 0, 
               min = 1, max = 1000, step = 1, width = NULL),
  numericInput(inputId = "n_groups", label = "How many groups for comparison does your data contain?", value = 0, 
  min = 1, max = 10, step = 1, width = NULL),
  numericInput(inputId = "min_fpkm", label = "Input the minimum FPKM value a gene must have (in at least
               one sample) to be used in this analysis:", value = 1, 
               min = 0, max = 10, step = 1, width = NULL),
  checkboxInput(inputId = "status_ok", label = "require STATUS = OK for gene to be includd in analysis?", value = TRUE,
                width = NULL),

  sidebarPanel(
    selectInput("gene_exp", "Download cleaned dataset", choices = "gene_exp", "test"),
    downloadButton('downloadData', 'Download'),
  ),
  mainPanel(
    tableOutput(output$gene_exp)
  )
  
  )))
  

server <- function(input, output) {
  #increase size of files that can be uploaded
  options(shiny.maxRequestSize = 5000*1024^2),
  
  #Load cuff diff
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

input$gene_exp <- gene.exp
datasetInput <- reactive({
  switch(datasetInput$gene_exp, "gene_exp" = gene.exp, "test" = test)
})
output$table <- renderTable({datasetInput()})
output$downloadData <- downloaadHandler(
  file_name = function() {paste(input$gene_exp, '.csv', sep = '') },
  content = function(file) { write.csv(datasetInput(), file)}
)

#
output$n_after <- renderPrint({ nrow(gene.exp) })
output$n_before <- renderPrint({nrow(gene.exp.full)})



}



shinyApp(ui = ui, server = server)








