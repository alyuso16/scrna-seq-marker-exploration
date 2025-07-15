library(shiny)
library(shinythemes)
library(shinyjs)

source("R/preprocessing.R")
source("R/umap_plots.R")
source("R/validation.R")

shinyApp(
  ui = fluidPage(theme = shinytheme("cosmo"),
    navbarPage("scRNa-seq",
      tabPanel("Tab 1",
        useShinyjs(),
        sidebarPanel(
          div(
            fileInput("file", "File upload:", accept = c(".h5")),
            actionButton("run", "Run", icon = icon("play"))
          )
        ),
        mainPanel(
          plotOutput("umap")
        )
      )
    )
  ),
  server <- function(input, output) {
    options(shiny.maxRequestSize = 1000 * 1024^2)

    disable("run")

    observeEvent(input$file, {
      if (!is.null(input$file)) {
        enable("run")
      } else {
        disable("run")
      }
    })

    seurat_obj <- eventReactive(input$run, {
      req(input$file)
      run_preprocessing(input$file$datapath, filetype = get_filetype(input$file))
    })

    output$umap <- renderPlot({
      req(seurat_obj())
      plot_umap(seurat_obj = seurat_obj())
    })
  }
)
