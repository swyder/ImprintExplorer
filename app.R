# Shiny web app for visualization of reanalyzed Arabidopsis datasets
# 
# Stefan Wyder
# April 2018

# TODO
# minimal read coverage not implemented yet
# color only stat. siginificant MEGs and PEGs

library(shiny)
library(ggplot2)
library(plotly)
library(stringr)


# Define UI for application
ui <- fluidPage(

   # App title
   titlePanel("ImprintExplorer"),

   sidebarLayout(
     sidebarPanel(
       selectInput("dataset", "Data set",
                   c("Hsieh", "Gehring", "Pignatta")),
       h3("Filters"),
       sliderInput(
         "MatProp",
         "% Maternal Reads",
         min = 0,
         max = 100,
         value = c(0, 100),
         step = 5
       ),
       sliderInput(
         "minCoverage",
         "Minimal Read Coverage",
         min = 0,
         max = 100,
         value = 1,
         step = 5
       ),
       sliderInput("fdr", "False Discovery Rate (FDR)",
                   0, 1, c(0, 1), step = 0.05),
       textInput("geneID", "Gene (e.g. AT4G00220)"
       )
      ),
     
      # Show plot and table tabs
      mainPanel(
        tabsetPanel(type = "tabs",
                    tabPanel("Plot",
                             h3("Reanalysis of Genomic Imprinting Studies in Arabidopsis"),
                             tags$a(href="https://doi.org/10.1101/180745", "Link to publication (bioRxiv 180745)"),
                             
                             plotlyOutput("scatterPlot"),     # plotOutput if not using plotly
                             
                             wellPanel(
                               span("Number of genes selected:",
                                    textOutput("n_genes")
                               ),
                               span("Number of MEGs (FDR<5%):",
                                    textOutput("n_MEGs_fdr5")
                             ),
                             span("Number of PEGs (FDR<5%):",
                                  textOutput("n_PEGs_fdr5")
                             ))),
                    tabPanel("ECDF",
                             h3("Empirical Cumulative Distribution"),
                             p("use without filtering"),
                             plotOutput("ECDF_Plot")
                             ),
                    tabPanel("Table", DT::dataTableOutput("table"),
                             br(),
                             downloadButton("downloadData", "Download csv table"))
        
      )
      )
   )
)



# Define server logic
server <- function(input, output) {
  # Supplementary table SI_Table_S1.csv from publication 
  allStudies <- read.table("SI_Table_S1.csv", header = TRUE, sep = "\t", skip = 2)
  colnames(allStudies)[seq(2,19)] <- paste(rep(c("Gehring", "Hsieh", "Pignatta"), each = 6), rep(c("logFC", "pval", "FDR", "type", "PercMaternal", "rank"), times = 3), sep = "_")
  colnames(allStudies)[seq(24,29)] <- paste("Wolff", c("logFC", "pval", "FDR", "type", "PercMaternal", "rank"), sep = "_")
  
  # create reactive data.frame as input to ggplotly
  dataInput <- reactive({
    # Select DataSet for filtering & plotting
    x <- switch(input$dataset,
                Pignatta = allStudies[, grepl("^Pignatta", colnames(allStudies))],
                Hsieh = allStudies[, grepl("^Hsieh", colnames(allStudies))],
                Gehring = allStudies[, grepl("^Gehring", colnames(allStudies))] 
    )
    
  # add common columns
    colnames(x) <- c("logFC", "pval", "FDR", "type", "PercMaternal", "rank")
    x <- cbind(x, allStudies[, c(1, 20, 21, 23)])
    colnames(x)[colnames(x) == 'Gene'] <- 'genes'
    
    # Apply filters
    m <- x %>%
      filter(
        logFC >= log2(input$MatProp[1]),
        logFC <= log2(input$MatProp[2]),
        FDR >= input$fdr[1],
        FDR <= input$fdr[2]
      )
    
    # Optional: filter by geneID
    if (!is.null(input$geneID) && input$geneID != "") {
      m <- m %>% filter(str_detect(input$geneID, as.character(genes)))
    }
    
    m <- as.data.frame(m)
    return(m)
  })
  
  # interactive plotting using plotly!
   output$scatterPlot <- renderPlotly({
     x <- dataInput() 
      p <- ggplot(data = x, aes(x = PercMaternal, y = FDR, text = genes, color = type, alpha = 0.5)) +   #, fill = clarity)) +
       geom_point() + xlab("% maternal reads") + ylab("FDR")
      ggplotly(p, tooltip = "genes") # to show multiple tooltip variables: tooltip=c("x", "y"))
   })
   
   # print number of filtered genes below plot
   output$n_genes <- renderText({ nrow(dataInput()) })
   output$n_MEGs_fdr5 <- renderText({ 
     dataInput() %>% filter(
       type == "MEG",
       FDR <= 0.05
     ) %>% nrow
     })
   output$n_PEGs_fdr5 <- renderText({
     dataInput() %>% filter(
       type == "PEG",
       FDR <= 0.05
     ) %>% nrow
   })
   
   output$ECDF_Plot <- renderPlot({
     ggplot(data=dataInput(), aes(PercMaternal)) + stat_ecdf(geom = "step", pad = FALSE) + xlab("% maternal read") + ylab("Proportion <= x")
   })

   output$table <- DT::renderDataTable(DT::datatable({ dataInput() }))
   
   # Downloadable csv of selected dataset ----
   output$downloadData <- downloadHandler(
     filename = function() {
       paste0(input$dataset, "_shiny.csv")
     },
     content = function(file) {
       write.csv(dataInput(), file, row.names = FALSE, sep = "\t")
     }
   )
}

# Run the application 
shinyApp(ui = ui, server = server)