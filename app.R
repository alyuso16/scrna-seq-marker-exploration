library(shiny)
library(shinythemes)
library(shinyjs)

source("R/preprocessing.R")
source("R/plots.R")
source("R/validation.R")

shinyApp(
  ui = fluidPage(theme = shinytheme("cosmo"),
    navbarPage("scRNa-seq",
      tabPanel("Tab 1",
        useShinyjs(),
        sidebarPanel(
          div(
            fileInput("file", "File upload:", multiple = FALSE, accept = c(".h5, .rds")),
            actionButton("run", "Run", icon = icon("play"))
          )
        ),
        mainPanel(
          plotOutput("plot"),
          uiOutput("plot_select"),
          uiOutput("marker_select")
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
      filetype <- get_filetype(input$file)

      if (filetype == "rds") {
        readRDS(input$file$datapath)
      } else {
        run_preprocessing(input$file)
      }
    })

    output$plot <- renderPlot({
      req(seurat_obj(), input$selected_plot)

      if (input$selected_plot == "UMAP") {
        plot_umap(seurat_obj())
      } else if (input$selected_plot == "Features") {
        req(input$selected_marker)
        plot_feature(seurat_obj(), feature = input$selected_marker)
      }
    })

    output$marker_select <- renderUI({
      req(seurat_obj())
      selectizeInput(
        "selected_marker",
        "Select marker gene:",
        choices = get_features(seurat_obj()),
        selected = get_features(seurat_obj())[1]
      )
    })

    output$plot_select <- renderUI({
      req(seurat_obj())
      selectizeInput(
        "selected_plot",
        "Select plot to display:",
        choices = c("UMAP", "Features"),
        selected = "UMAP"
      )
    })
  }
)
