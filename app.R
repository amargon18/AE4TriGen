library(shiny)
library(bslib)
library(shinybusy)
library(shinyjs)
library(shinyFiles)
library(plotly)
library(httr)
library(jsonlite)
library(DT)
source("helpers.R")

theme <- bs_theme(
  version    = 5,
  bootswatch = "lux",
  primary    = "#DC143C",
  secondary  = "#8B0000"
)

ui <- fluidPage(
  useShinyjs(),
  theme = theme,
  add_busy_spinner(
    spin     = "fading-circle",
    position = "bottom-right",
    margins  = c(20,50),
    color    = "#DC143C"
  ),
  tags$head(
    tags$link(rel = "icon", type = "image/png", href = "logoAE4TriGen.png"),
    tags$script(HTML("
      document.addEventListener('DOMContentLoaded', function() {
        var t = [].slice.call(
          document.querySelectorAll('[data-bs-toggle=\"tooltip\"]')
        );
        t.forEach(el => new bootstrap.Tooltip(el));
      });
      .container-fluid {
        padding-left: 0;
        padding-right: 0;
      };
    "))
  ),
  
fluidRow(
  style = "background: linear-gradient(90deg, #DC143C, #8B0000);",
  column(8,
         div(
           style = "text-align:center; margin:20px 0;padding: 15px; border-radius: 10px;",
           tags$h1(
             "ARRAYEXPRESS PROCESSOR for TriGen Algorithm",
             style = "
              font-size:1.2rem; 
              font-weight:bold;
              white-space: nowrap;  /* Evita el salto de línea */
              overflow: hidden;  /* Previene desbordamiento */
              color: white;  /* Letras en blanco */
              font-family: 'Montserrat', sans-serif;
            "
           )
         )
  ),
  column(4,
         div(
           style = "text-align:right; margin:20px 0;",
           actionButton(
             inputId = "open_search",
             label   = NULL,
             icon    = icon("search"),
             class   = "btn btn-link search-icon"
           )
         )
  )
),

# Estilos adicionales para mejorar el diseño
tags$style("
  .search-icon:hover {
    transform: scale(1.1);
    transition: 0.3s ease-in-out;
  }

  .btn-link:hover {
    color: #DC143C;
    transition: 0.3s;
  }

  
 
  
"),


  hr(),
  
  uiOutput("search_table"),
  
  # SELECTOR DE DIRECTORIO
  fluidRow(
    column(6,
           shinyDirButton(
             "choose_dir", "Seleccionar directorio",
             "Selecciona carpeta base", FALSE
           )
    ),
    column(6,
           verbatimTextOutput("chosen_dir")
    )
  ),
  hr(),
  
  # MAIN LAYOUT
  fluidRow(
    column(4,
           wellPanel(
             selectizeInput("accession", "Experimento:", choices=NULL,
                            options=list(create=TRUE), width="100%"),
             actionButton("download_btn","Descargar datos",
                          icon=icon("download"),class="btn btn-danger w-100 mb-2"),
             actionButton("process_btn","Crear archivos",
                          icon=icon("cogs"),class="btn btn-danger w-100 mb-2"),
             actionButton("runtrigen_btn","Ejecutar TriGen",
                          icon=icon("play"),class="btn btn-danger w-100 mb-2"),
             selectInput("genome","Genoma (PANTHER):",choices=NULL),
             selectInput("annot_dataset","Dataset anotación:",choices=NULL),
             actionButton("enrich_btn","Enriquecimiento",
                          icon=icon("chart-bar"),class="btn btn-danger w-100 mb-2"),
             actionButton("consultar_genes","Info genes",
                          icon=icon("dna"),class="btn btn-danger w-100 mb-2"),
             actionButton("expr_btn","Expr. interactiva",
                          icon=icon("chart-line"),class="btn btn-danger w-100 mb-2"),
             actionButton("metrics_btn","Métricas",
                          icon=icon("table"),class="btn btn-danger w-100")
           )
    ),
    column(8,
           navset_pill( 
             tabPanel("Estado",
                      verbatimTextOutput("download_status"),
                      verbatimTextOutput("process_status")
             ),
             tabPanel("Archivos",
                      verbatimTextOutput("archivos_creados")
             ),
             tabPanel("Enriquecimiento",
                      uiOutput("plots_enrichment")
             ),
             tabPanel("Genes",
                      uiOutput("geneinfo_ui")
             ),
             tabPanel("Expresión",
                      uiOutput("expr_plots")
             ),
             tabPanel("Métricas",
                      dataTableOutput("metrics_table")
             )
           )
    )
  ),
  fluidRow(
    class = "footer bg-light border-top",
    style = "padding: 1rem 0; font-size: 0.85rem;",
    
    # Autor
    column(
      4,
      tags$ul(
        class = "list-unstyled mb-0",
        tags$li(
          class = "text-muted",
          "Desarrollado por Alejandro Márquez González"
        )
      )
    ),
    
    # Enlace
    column(
      4,
      div(
        class = "text-center",
        tags$a(
          href   = "https://www.ebi.ac.uk/arrayexpress/",
          target = "_blank",
          class  = "text-primary fw-bold text-decoration-none",
          "Visita ArrayExpress"
        )
      )
    ),
    
    # Logos
    column(
      4,
      div(
        class = "d-flex justify-content-end align-items-center",
        img(src = "logo2.png", height = "40px")
      )
    )
  )
  
)

server <- function(input, output, session) {
  roots <- c(Home = normalizePath("~"))
  
  shinyDirChoose(
    input,
    id    = "choose_dir",
    roots = roots,
    session = session,
    filetypes = c("", "txt", "csv")
  )
  
  baseDir <- reactive({
    if (is.null(input$choose_dir) || length(input$choose_dir) == 0) return(NULL)
    normalizePath(
      parseDirPath(roots, input$choose_dir),
      winslash = "/")
  })
  
  output$chosen_dir <- renderText({
    bd <- baseDir()
    if (!shiny::isTruthy(bd)) {
      "Debes seleccionar el directorio"
    } else {
      paste("Directorio seleccionado:", bd)
    }
  })
  
  # 1) Poblar accesiones
  observe({
    acc <- list.dirs(file.path(baseDir(),"DatosExp"),
                     full.names=FALSE, recursive=FALSE)
    acc <- acc[nzchar(acc) & acc!="."]
    updateSelectInput(session, "accession", choices=acc, selected=acc[1])
  })
  
  # 2) PANTHER dropdowns
  observe({
    resp <- POST("http://pantherdb.org/services/oai/pantherdb/supportedgenomes",
                 accept_json())
    dat  <- fromJSON(httr::content(resp, "text", encoding="UTF-8"))
    gen  <- dat$search$output$genomes$genome$taxon_id
    updateSelectInput(session, "genome", choices=gen, selected=gen[1])
  })
  observe({
    resp <- GET("http://pantherdb.org/services/oai/pantherdb/supportedannotdatasets",
                accept_json())
    dat  <- fromJSON(httr::content(resp, "text", encoding="UTF-8"))
    ds   <- dat$search$annotation_data_sets$annotation_data_type$id
    updateSelectInput(session, "annot_dataset", choices=ds, selected=ds[1])
  })
  
  # 3) Descargar datos
  observeEvent(input$download_btn, {
    req(baseDir(), input$accession)
    output$download_status <- renderText("Comprobando y descargando…")
    result <- downloadData(input$accession, base_dir = baseDir())
    output$download_status <- renderText({
      if (identical(result, "already_downloaded")) {
        files <- list.files(file.path(baseDir(),"DatosExp", input$accession),
                            full.names=TRUE, recursive=TRUE)
        paste0("Este experimento fue descargado previamente. Se encontraron los siguientes archivos:\n", 
               paste(files, collapse="\n"))
      } else if (is.null(result)) {
        paste("Error al descargar:", input$accession)
      } else {
        files <- list.files(file.path(baseDir(),"DatosExp", input$accession),
                            full.names=TRUE, recursive=TRUE)
        paste0("Descarga completa:\n", paste(files, collapse="\n"))
      }
    })
  })
  
  # 4) Crear archivos
  observeEvent(input$process_btn, {
    req(baseDir(), input$accession)
    archivos <- tryCatch(fileCreationFinal(input$accession,
                                           base_dir = baseDir()), 
                         error=function(e) NULL)
    archivos <- list.files(file.path(baseDir(),"TrLab5/Resources", input$accession),
                           full.names=TRUE, recursive=TRUE)
    output$archivos_creados <- renderText(paste(archivos, collapse="\n"))
  })
  
  # 5) Ejecutar TriGen
  observeEvent(input$runtrigen_btn, {
    req(baseDir(),input$accession)
    output$process_status <- renderText("Ejecutando TriGen…")
    out <- tryCatch(runTriGen(input$accession,
                              base_dir = baseDir()), error=function(e) NULL)
    output$process_status <- renderText(
      if (is.null(out)) "Error al ejecutar TriGen"
      else paste("TriGen ejecutado. Resultados en:", out$dir)
    )
  })
  
  # 6) Enriquecimiento funcional
  observeEvent(input$enrich_btn, {
    req(baseDir(),input$accession, input$genome, input$annot_dataset)
    show_modal_spinner(text="Analizando…", color="#DC143C")
    genes_dir <- file.path(baseDir(), "TrLab5/OutputTriGen",
                           input$accession, input$accession, "genes")
    txts <- list.files(genes_dir, pattern="\\.txt$", full.names=TRUE)
    if (!length(txts)) {
      remove_modal_spinner()
      output$plots_enrichment <- renderUI(div(class="alert alert-warning",
                                              "No se encontraron archivos de genes."))
      return()
    }
    plots <- lapply(txts, function(f) {
      name <- tools::file_path_sans_ext(basename(f))
      genes <- readLines(f, warn=FALSE)
      df    <- runOverrepAnalysis(genes, input$genome, input$annot_dataset)
      top   <- getTopTerms(df, n=10)
      plotEnrichment(top, title=name)
    })
    remove_modal_spinner()
    output$plots_enrichment <- renderUI({
      tabs <- lapply(seq_along(plots), function(i) {
        nm <- tools::file_path_sans_ext(basename(txts[i]))
        tabPanel(nm, plotlyOutput(paste0("plt_", i)))
      })
      do.call(tabsetPanel, tabs)
    })
    for (i in seq_along(plots)) {
      local({
        idx <- i
        output[[paste0("plt_", idx)]] <- renderPlotly(plots[[idx]])
      })
    }
  })
  
  # 7) Consulta de genes
  observeEvent(input$consultar_genes, {
    req(baseDir(), input$accession, input$genome)
    
    res_list <- queryGeneInfo(input$accession,
                              input$genome, base_dir = baseDir())
    
  
    if (all(sapply(res_list, function(x) is.null(x) || nrow(x) == 0))) {
      output$geneinfo_ui <- renderUI({
        div(class="alert alert-warning",
            "No se encontraron resultados válidos.")
      })
      return()
    }
    
    # 1) Generamos las pestañas
    output$geneinfo_ui <- renderUI({
      tabs <- lapply(names(res_list), function(name) {
        tabPanel(
          title = name,
          dataTableOutput(outputId = paste0("geneinfo_tbl_", name))
        )
      })
      do.call(tabsetPanel, tabs)
    })
    
    # 2) Para cada data.frame creamos su renderDataTable
    for (name in names(res_list)) {
      local({
        nm <- name
        df <- res_list[[nm]]
        out_id <- paste0("geneinfo_tbl_", nm)
        output[[out_id]] <- renderDataTable({
          df
        }, options = list(pageLength = 10, scrollX = TRUE))
      })
    }
  })
  
  # 8) Expresión interactiva
  observeEvent(input$expr_btn, {
    req(baseDir(),input$accession)
    plots <- geneExpInteractive(input$accession, base_dir = baseDir())
    if (!length(plots)) {
      output$expr_plots <- renderUI(
        div(class="alert alert-warning","No hay archivos de coordenadas.")
      )
      return()
    }
    output$expr_plots <- renderUI({
      top_tabs <- lapply(names(plots), function(fn) {
        inner <- lapply(names(plots[[fn]]), function(cond) {
          pid <- paste0("expr_", fn, "_", cond)
          tabPanel(cond, plotlyOutput(pid))
        })
        tabPanel(fn, do.call(tabsetPanel, inner))
      })
      do.call(tabsetPanel, top_tabs)
    })
    for (fn in names(plots)) {
      for (cond in names(plots[[fn]])) {
        local({
          pod <- paste0("expr_", fn, "_", cond)
          pobj <- plots[[fn]][[cond]]
          output[[pod]] <- renderPlotly(pobj)
        })
      }
    }
  })
  
  # 9) Métricas TriGen
  observeEvent(input$metrics_btn, {
    req(baseDir(), input$accession)
    df <- showMetrics(input$accession, base_dir = baseDir())
    output$metrics_table <- renderDataTable(df,
                                            options=list(pageLength=10,
                                                         scrollX=TRUE))
  })
  
  observeEvent(input$open_search, {
    showModal(modalDialog(
      title = "Buscar experimentos (BioStudies)",
      textInput("search_term", "Experimental design:", value = ""),
      actionButton("search_btn", "Buscar", icon = icon("search")),
      hr(),
      # Aquí es donde la tabla se va a dibujar
      tableOutput("search_table"),
      footer = modalButton("Cerrar")
    ))
  })
  
  
  observeEvent(input$search_btn, {
    req(input$search_term)  
    # Llama a la función que definiste en helpers.R:
    resultados <- consultaExp(input$search_term)
    
    if (nrow(resultados) == 0) {
      # Si no hay filas, muestra un mensaje “vacío”
      output$search_table <- renderTable({
        data.frame(Mensaje = "No se encontraron experimentos.")
      }, rownames = FALSE, colnames = FALSE)
    } else {
      # Si sí hay filas, renderiza las columnas que te interesen
      output$search_table <- renderTable({
        # Por ejemplo, mostramos solo “accession” y “title”
        resultados[, c("accession", "title", "files"), drop = FALSE]
      }, rownames = FALSE,
      hover = TRUE,
      border = TRUE,)
    }
  })
  
  
}


shinyApp(ui, server)
