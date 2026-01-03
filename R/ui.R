library(shinythemes)
library(shinyjs)
library(bslib)
library(DT)

ui <- page_navbar(
  title = "scRNA-seq Marker Explorer",
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
        dataTableOutput("marker_table")
      )
    ),
    mainPanel(
      fluidRow(
        column(
          width = 9,
          plotOutput("umap_plot", height = "600px")
        )
      ),
      fluidRow(
        column(
          width = 9,
          plotOutput("feature_plot", height = "600px")
        )
      ),
      fluidRow(
        column(
          width = 9,
          plotOutput("marker_violin_plot", height = "600px")
        )
      )
    )
  )
)