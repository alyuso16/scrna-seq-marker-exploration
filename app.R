library(shiny)
library(shinythemes)

source("R/preprocessing.R")
source("R/umap_plots.R")

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
      filetype <- tolower(tools::file_ext(input$file$name))
      validate(need(filetype %in% c("h5"), "Upload .h5 file"))

      withProgress(message = "Running preprocessing pipeline...", value = 0, {
        progress <- shiny::Progress$new()
        on.exit(progress$close())

        run_preprocessing(input$file$datapath, filetype = filetype, progress = progress)
      })
    })

    output$umap <- renderPlot({
      req(seurat_obj())
      plot_umap(seurat_obj = seurat_obj())
    })
  }
)
