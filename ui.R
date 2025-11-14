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
          fileInput("file", "File upload:", multiple = FALSE, accept = c(".h5", ".h5ad", ".rds")),
          uiOutput("upload_msg"),
        ),
        div(
          style = "margin-bottom: 20px;",
          actionButton("run", "Run", icon = icon("play")),
        ),

        textInput("object_save_name", "Save Seurat object as:"),
        downloadButton("download_object", "Download Seurat Object"),
      ),
      mainPanel(
        fluidRow(
          column(
            width = 9,
            plotOutput("umap_plot")
          ),
          column(
            width = 3,
            downloadButton("download_umap", "Download UMAP plot")
          )
        ),
        fluidRow(
          column(
            width = 9,
            plotOutput(("feature_plot"))
          ),
          column(
            width = 3,
            uiOutput("marker_select"),
            downloadButton("download_feature_plot", "Download this feature plot")
          )
        ),
        fluidRow(
          column(
            width = 9,
            plotOutput(("marker_violin_plot"))
          ),
          column(
            width = 3,
            downloadButton("download_violin_plot", "Download this violin plot")
          )
        )
      )
    )
  )
)