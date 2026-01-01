library(shinythemes)
library(shinyjs)
library(bslib)

ui <- page_navbar(
  title = "scRNA-seq",
  useShinyjs(),
  sidebarLayout(
    sidebarPanel(
      div(
        style = "margin-bottom: 20px;",
        fileInput("file", "RDS file upload:", multiple = FALSE, accept = c(".rds"))
      ),
      div(
        style = "margin-bottom: 20px;",
        input_task_button("run", "Run", icon = icon("play")),
      ),
      div(
        style = "margin-bottom: 20px;",
        uiOutput("cluster_select"),
        uiOutput("marker_sort"),
        uiOutput("marker_select")
      )
    ),
    mainPanel(
      fluidRow(
        column(
          width = 9,
          plotOutput("umap_plot")
        ),
        column(
          width = 3,
          uiOutput("download_umap_button")
        )
      ),
      fluidRow(
        column(
          width = 9,
          plotOutput(("feature_plot"))
        ),
        column(
          width = 3,
          uiOutput("download_feature_button")
        )
      ),
      fluidRow(
        column(
          width = 9,
          plotOutput(("marker_violin_plot"))
        ),
        column(
          width = 3,
          uiOutput("download_violin_button")
        )
      )
    )
  )
)