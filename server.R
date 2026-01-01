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
    req(seurat_obj(), input$selected_marker)
    plot_feature(seurat_obj(), marker = input$selected_marker)
  })

  marker_violin_plot <- reactive({
    req(seurat_obj(), input$selected_marker)
    plot_marker_violin(seurat_obj(), marker = input$selected_marker)
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

  output$cluster_select <- renderUI({
    req(seurat_obj())
    selectizeInput(
      "selected_cluster",
      "Filter markers by cluster:",
      choices = c(list("All Clusters"), get_cluster_labels(seurat_obj())),
      selected = "All Clusters"
    )
  })

  output$marker_sort <- renderUI({
    req(seurat_obj())
    selectizeInput(
      "selected_sort",
      "Sort markers:",
      choices = list(
        "Expression: avg_log2FC" = 1,
        "Expression: pct.1" = 2
      ),
      selected = 1
    )
  })

  output$marker_select <- renderUI({
    req(seurat_obj(), all_markers(), input$selected_cluster, input$selected_sort)

    selectizeInput(
      "selected_marker",
      "Select marker gene:",
      choices = filter_markers(
        all_markers(),
        input$selected_cluster,
        input$selected_sort
      ),
      selected = filter_markers(
        all_markers(),
        input$selected_cluster,
        input$selected_sort
      )[1]
    )
  })

  output$download_umap_button <- renderUI({
    req(umap_plot())
    downloadButton("download_umap", "Download UMAP plot")
  })

  output$download_umap <- downloadHandler(
    filename = function() {
      paste0(trimws(input$file$name), "_umap.png")
    },
    content = function(file) {
      ggsave(file, plot = umap_plot(), device = "png", width = 20, height = 15)
    }
  )

  output$download_feature_plot <- downloadHandler(
    filename = function() {
      paste0(trimws(input$file$name), "_feature_plot_", input$selected_marker, ".png")
    },
    content = function(file) {
      ggsave(file, plot = feature_plot(), device = "png", width = 20, height = 15)
    }
  )

  output$download_violin_button <- renderUI({
    req(marker_violin_plot())
    downloadButton("download_violin_plot", "Download violin plot")
  })

  output$download_violin_plot <- downloadHandler(
    filename = function() {
      paste0(trimws(input$file$name), "_violin_plot_", input$selected_marker, ".png")
    },
    content = function(file) {
      ggsave(file, plot = marker_violin_plot(), device = "png")
    }
  )
}