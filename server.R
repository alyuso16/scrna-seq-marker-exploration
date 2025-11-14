library(shiny)
library(shinythemes)
library(shinyjs)
library(bslib)

source("R/preprocessing.R")
source("R/plots.R")
source("R/validation.R")

server <- function(input, output) {
  options(shiny.maxRequestSize = 2000 * 1024^2)

  disable("run")
  disable("download_umap")
  disable("download_feature_plot")

  observeEvent(input$file, {
    if (!is.null(input$file)) {
      enable("run")
    } else {
      disable("run")
    }
  })

  # observeEvent(input$file, {
  #   if (!is.null(input$file)) {
  #     enable("run")
  #   } else {
  #     disable("run")
  #   }
  # })

  # observeEvent(input$file, {
  #   if (!is.null(input$file)) {
  #     enable("run")
  #   } else {
  #     disable("run")
  #   }
  # })

  seurat_obj <- eventReactive(input$run, {
    req(input$file)
    filetype <- get_filetype(input$file)

    if (filetype == "rds") {
      readRDS(input$file$datapath)
    } else {
      run_preprocessing(input$file)
    }
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
    req(seurat_obj())
    plot_marker_violin(seurat_obj(), marker = input$selected_marker)
  })

  output$upload_msg <- renderUI({
    req(input$file)
    filetype <- get_filetype(input$file)

    tags$span(
      style = "color: green;",
      paste(
        "Uploaded ", filetype, " file,",
        if (filetype == "rds") {
          "run to skip preprocessing"
        } else if (filetype == "h5" | filetype == "h5ad") {
          "run to begin preprocessing"
        }
      )
    )

  })

  output$download_object <- downloadHandler(
    filename = function() {
      paste0(trimws(input$object_save_name), ".rds")
    },
    content = function(file) {
      req(seurat_obj())
      saveRDS(seurat_obj(), file)
    }
  )

  output$umap_plot <- renderPlot({
    umap_plot()
  })

  output$feature_plot <- renderPlot({
    feature_plot()
  })

  output$marker_violin_plot <- renderPlot({
    marker_violin_plot()
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

  output$download_umap <- downloadHandler(
    filename = function() {
      paste0(trimws(input$file$name), "_umap.png")
    },
    content = function(file) {
      req(umap_plot())
      ggsave(file, plot = umap_plot(), device = "png", width = 20, height = 15)
    }
  )

  output$download_feature_plot <- downloadHandler(
    filename = function() {
      paste0(trimws(input$file$name), "_feature_plot_", input$selected_marker, ".png")
    },
    content = function(file) {
      req(feature_plot())
      ggsave(file, plot = feature_plot(), device = "png", width = 20, height = 15)
    }
  )

  output$download_violin_plot <- downloadHandler(
    filename = function() {
      paste0(trimws(input$file$name), "_violin_plot_", input$selected_marker, ".png")
    },
    content = function(file) {
      req(marker_violin_plot())
      ggsave(file, plot = marker_violin_plot(), device = "png")
    }
  )
}