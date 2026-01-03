library(shinythemes)
library(shinyjs)
library(bslib)

source("R/plots.R")
source("R/helper.R")

server <- function(input, output) {
  options(shiny.maxRequestSize = 2000 * 1024^2)

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
    }
  })

  all_markers <- reactive({
    req(seurat_obj())

    withProgress(
      message = "Calculating clusters...",
      value = 0,
      {
        incProgress(0.5)
        features <- get_features(seurat_obj())
        incProgress(0.7)
        features
      }
    )
  })

  umap_plot <- reactive({
    req(seurat_obj())
    plot_umap(seurat_obj())
  })

  feature_plot <- reactive({
    req(seurat_obj())

    if (is.null(input$marker_table_cell_clicked$value)) {
      plot_feature(seurat_obj(), marker = all_markers()$gene[1])
    } else {
      plot_feature(seurat_obj(), marker = all_markers()$gene[input$marker_table_cell_clicked$row])
    }
  })

  marker_violin_plot <- reactive({
    req(seurat_obj())

    if (is.null(input$marker_table_cell_clicked$value)) {
      plot_marker_violin(seurat_obj(), marker = all_markers()$gene[1])
    } else {
      plot_marker_violin(seurat_obj(), marker = all_markers()$gene[input$marker_table_cell_clicked$row])
    }
  })

  output$umap_plot <- renderPlot({
    umap_plot()
  })

  output$feature_plot <- renderPlot({
    feature_plot()
  })

  output$marker_violin_plot <- renderPlot({
    marker_violin_plot()
  })

  output$marker_table <- renderDataTable({
    req(seurat_obj(), all_markers())

    datatable(
      all_markers(),
      selection = "single",
      filter = "top"
    )
  })
}