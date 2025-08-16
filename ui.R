library(shiny)
library(shinythemes)
library(shinyjs)
library(bslib)

ui <- page_navbar(
  title = "scRNA-seq",
  useShinyjs(),
  tabPanel("Tab 1",
    sidebarLayout(
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
)