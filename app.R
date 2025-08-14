library(shiny)
library(shinythemes)
library(shinyjs)
library(bslib)

source("R/preprocessing.R")
source("R/plots.R")
source("R/validation.R")

ui = page_navbar(
  title = "scRNA-seq",
  useShinyjs(),
  tabPanel("Tab 1",
    sidebarPanel(
      div(
        style = "margin-bottom: 20px;",
        fileInput("file", "File upload:", multiple = FALSE, accept = c(".h5, .rds")),
        uiOutput("upload_msg"),
      ),

      div(
        style = "margin-bottom: 20px;",
        actionButton("run", "Run", icon = icon("play")),
      ),

      textInput("save_name", "Save Seurat object as:"),
      downloadButton("download_object", "Download Seurat Object"),
    ),
    mainPanel(
      fluidRow(
        column(
          width = 9,
          plotOutput("umap_plot")
        ),
      ),
      fluidRow(
        column(
          width = 9,
          plotOutput(("feature_plot"))
        ),
        
        column(
          width = 3,
          uiOutput("marker_select")
        )
      )
    )
  )
)

# ui = fluidPage(
#   navbarPage("scRNA-seq",
#     tabPanel("Tab 1",
#       useShinyjs(),
#       sidebarPanel(
#         div(
#           style = "margin-bottom: 20px;",
#           fileInput("file", "File upload:", multiple = FALSE, accept = c(".h5, .rds")),
#           uiOutput("upload_msg"),
#         ),
#         div(
#           style = "margin-bottom: 20px;",
#           actionButton("run", "Run", icon = icon("play")),
#         ),

#         textInput("save_name", "Save Seurat object as:"),
#         downloadButton("download_object", "Download Seurat Object"),
#       ),
#       mainPanel(
#         fluidRow(
#           column(
#             width = 9,
#             plotOutput("umap_plot")
#           ),
#         ),
#         fluidRow(
#           column(
#             width = 9,
#             plotOutput(("feature_plot"))
#           ),
          
#           column(
#             width = 3,
#             uiOutput("marker_select")
#           )
#         )
#       )
#     )
#   )
# )

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

  output$upload_msg <- renderUI({
    req(input$file)
    filetype <- get_filetype(input$file)

    tags$span(
      style = "color: green;",
      paste(
        "Uploaded ", filetype, " file,",
        if (filetype == "rds") {
          "run to skip preprocessing..."
        } else if (filetype == "h5") {
          "run to begin preprocessing..."
        }
      )
    )
    
  })

  output$download_object <- downloadHandler(
    filename = function() {
      paste0(trimws(input$save_name), ".rds")
    },
    content = function(file) {
      req(seurat_obj())
      saveRDS(seurat_obj(), file)
    }
  )

  output$umap_plot <- renderPlot({
    req(seurat_obj())

    plot_umap(seurat_obj())
  })

  output$feature_plot <- renderPlot({
    req(seurat_obj(), input$selected_marker)

    plot_feature(seurat_obj(), marker = input$selected_marker)
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
}

shinyApp(ui, server)