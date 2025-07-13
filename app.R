library(shiny)
library(shinythemes)

source("R/preprocessing.R")
source("R/umap_plots.R")
source("R/validation.R")

options(shiny.maxRequestSize = 1000 * 1024^2)

 
shinyApp(
  ui = fluidPage(theme = shinytheme("cosmo"),
    navbarPage("scRNa-seq",
      tabPanel("Tab 1",
        sidebarPanel(
          fileInput("file", "File upload:", accept = c(".h5"))
        ),
        mainPanel(
          plotOutput("umap")
        )
      )
    )
  ),
  server <- function(input, output) {
    seurat_obj <- reactive({
      req(input$file)
      run_preprocessing(input$file$datapath, filetype = get_filetype(input$file))
    })

    output$umap <- renderPlot({
      req(seurat_obj())
      plot_umap(seurat_obj = seurat_obj())
    })
  }
)
