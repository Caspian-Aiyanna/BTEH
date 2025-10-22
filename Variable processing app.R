library(shiny)
library(shinydashboard)
library(sf)
library(raster)
library(terra)
library(DT)
library(leaflet)
library(shinyFiles)

ui <- dashboardPage(
  dashboardHeader(title = "Raster Processing for 
                  Uniform Resolution, Extent and CRS"),
  dashboardSidebar(collapsed = TRUE),
  dashboardBody(
    fluidRow(
      box(
        title = "Input Files",
        width = 16,
        fileInput("shp_files", 
                  "Upload Shapefile (select all .shp, .shx, .dbf, .prj files)", 
                  multiple = TRUE,
                  accept = c('.shp','.shx','.dbf','.prj')),
        
        fileInput("ref_raster", 
                  "Upload Reference Raster (.tif)", 
                  multiple = FALSE,
                  accept = ".tif"),
        
        fileInput("raster_files", 
                  "Upload Raster Files to Process (.tif)", 
                  multiple = TRUE,
                  accept = ".tif"),
        
        shinyDirButton("dir", "Choose Output Directory", "Please select a directory"),
        verbatimTextOutput("dir_path"),
        
        selectInput("method", 
                    "Resampling Method:",
                    choices = c("bilinear", "nearest", "cubic")),
        
        actionButton("process_btn", 
                     "Process Rasters", 
                     class = "btn-primary"),
        hr(),
        DTOutput("results_table")
      )
    )
  )
)

server <- function(input, output, session) {
  options(shiny.maxRequestSize = 0)
  
  volumes <- getVolumes()()
  shinyDirChoose(input, "dir", roots = volumes, session = session)
  
  dir_path <- reactive({
    if (!is.null(input$dir)) {
      parseDirPath(volumes, input$dir)
    }
  })
  
  output$dir_path <- renderText({
    if (!is.null(dir_path())) {
      dir_path()
    }
  })
  
  rv <- reactiveValues(
    shapefile = NULL,
    reference_raster = NULL,
    results = NULL
  )
  
  observeEvent(input$shp_files, {
    req(input$shp_files)
    
    temp_dir <- tempdir()
    sapply(input$shp_files$datapath, function(x) {
      file.copy(x, file.path(temp_dir, basename(input$shp_files$name[which(input$shp_files$datapath == x)])))
    })
    
    shp_file <- list.files(temp_dir, pattern = "\\.shp$", full.names = TRUE)[1]
    if (!is.null(shp_file)) {
      rv$shapefile <- st_read(shp_file, quiet = TRUE)
    }
  })
  
  observeEvent(input$ref_raster, {
    req(input$ref_raster)
    rv$reference_raster <- terra::rast(input$ref_raster$datapath)
  })
  
  observeEvent(input$process_btn, {
    req(rv$shapefile, rv$reference_raster, input$raster_files, dir_path())
    
    dir.create(dir_path(), showWarnings = FALSE, recursive = TRUE)
    
    results_df <- data.frame(
      Filename = character(),
      Status = character(),
      Resolution = numeric(),
      stringsAsFactors = FALSE
    )
    
    withProgress(message = 'Processing rasters', value = 0, {
      for(i in 1:length(input$raster_files$datapath)) {
        tryCatch({
          input_rast <- terra::rast(input$raster_files$datapath[i])
          
          processed_rast <- resample(input_rast, rv$reference_raster, method = input$method)
          processed_rast <- mask(processed_rast, rv$shapefile)
          processed_rast <- crop(processed_rast, rv$shapefile)
          
          output_file <- file.path(dir_path(), 
                                   paste0("processed_", input$raster_files$name[i]))
          writeRaster(processed_rast, output_file, overwrite = TRUE)
          
          results_df <- rbind(results_df, 
                              data.frame(Filename = input$raster_files$name[i],
                                         Status = "Success",
                                         Resolution = res(rv$reference_raster)[1]))
          
        }, error = function(e) {
          results_df <- rbind(results_df,
                              data.frame(Filename = input$raster_files$name[i],
                                         Status = paste("Error:", e$message),
                                         Resolution = NA))
        })
        incProgress(1/length(input$raster_files$datapath))
      }
    })
    
    rv$results <- results_df
  })
  
  output$results_table <- renderDT({
    req(rv$results)
    rv$results
  })
}

shinyApp(ui, server)
