library(httr)
library(jsonlite)
library(dbplyr)
library(biomaRt)
library(ArrayExpress)
library(affy)
library(xml2)
library(RCurl)
library(XML)
library(tidyr)
library(ggplot2)
library(plotly)

### DESCARGA DE DATOS
downloadData <- function(accesion_number, base_dir) {
  # Construccion directorio local para almacenamiento
  datos_dir <- file.path(base_dir, "DatosExp")
  if (!dir.exists(datos_dir)) {
    dir.create(datos_dir, recursive = TRUE)
  }
  
  # Creacion del directorio en caso de no haberse creado previamente
  path_folder <- file.path(datos_dir, accesion_number)
  if (dir.exists(path_folder)) {
    archivos_existentes <- list.files(path_folder, recursive = TRUE)
    if (length(archivos_existentes) > 0) {
      message("El experimento ya ha sido descargado previamente.")
      return("already_downloaded")
    }
  } else {
    dir.create(path_folder, recursive = TRUE)
  }
  
  # Acceso al experimento y posterior descarga
  # Uso funcion auxiliar en caso de error.
  experiment <- tryCatch({
    getAE(accesion_number, type = "full", path = path_folder)
    path_folder 
  }, error = function(e) {
    message("Error en la descarga: ", e$message)
    downloadFromFTP(accesion_number, dest_dir = datos_dir)
  })
  
  return(experiment)
}

downloadFromFTP <- function(accession, dest_dir, unzip_files = TRUE) {
  ftp_base  <- "ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment"
  exp_type  <- toupper(sub("(^E-)([^-]+).*", "\\2", accession))
  ftp_url   <- file.path(ftp_base, exp_type, accession)
  local_dir <- file.path(dest_dir, accession)
  if (!dir.exists(local_dir)) dir.create(local_dir, recursive = TRUE)
  
  zip_filename   <- paste0(accession, ".raw.1.zip")
  zip_url        <- file.path(ftp_url, zip_filename)
  local_zip_path <- file.path(local_dir, zip_filename)
  
  tryCatch({
    download.file(zip_url, local_zip_path, mode = "wb")
    if (unzip_files) unzip(local_zip_path, exdir = local_dir)
  }, error = function(e) {
    message("No se pudo descargar ZIP: ", e$message)
  })
  
  for (ext in c(".sdrf.txt", ".idf.txt", ".adf.txt")) {
    fname <- paste0(accession, ext)
    url   <- file.path(ftp_url, fname)
    dest  <- file.path(local_dir, fname)
    tryCatch({
      download.file(url, dest, mode = "wb")
    }, error = function(e) {
      message("No encontrado ", fname)
    })
  }
  
  return(local_dir)
}

### 2. PROCESAMIENTO Y CREACIÓN DE ARCHIVOS

fileCreationFinal <- function(accesion_number, base_dir){
 
  ### FICHERO ARCHIVO PUNTOS DE TIEMPO
  sdrf_path <- file.path(base_dir,"DatosExp", accesion_number, paste0(accesion_number, ".sdrf.txt"))
  sdrf_data <- read.delim(sdrf_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  columnas_tiempo <- grep("time", colnames(sdrf_data), value = TRUE, ignore.case = TRUE)
  
  output_dir <- file.path(base_dir, "TrLab5/resources", accesion_number)
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  output_path <- file.path(output_dir, "time_points.txt")
  
  time_columns_int <- columnas_tiempo[sapply(sdrf_data[columnas_tiempo], is.integer)]
  
  sdrf_agrupado <- sdrf_data %>%
    dplyr::group_by(sdrf_data[[time_columns_int]]) %>% 
    dplyr::summarise(Scans = paste(Array.Data.File, collapse = ", ")) %>%
    dplyr::arrange(time_columns_int)  # Ordenar por tiempo
  
  chr_column <- sdrf_agrupado[[2]]
  counts <- sapply(strsplit(chr_column, ","), length)
  most_common_count <- as.numeric(names(sort(table(counts), decreasing = TRUE)[1]))
  valid_rows <- counts == most_common_count
  cleaned_sdrf_agrupado <- sdrf_agrupado[valid_rows, ]
  valores_tiempo <- cleaned_sdrf_agrupado[[1]]
  time_unit <- unique(sdrf_data[[columnas_tiempo[[2]]]])
  cleaned_value_times <- paste(valores_tiempo, time_unit)
  
  writeLines(cleaned_value_times, con = output_path)
  
  
  ### ARCHIVO CONDICIONES EXPERIMENTALES y MATRICES DE EXPRESIÓN POR PUNTO TEMPORAL.
  cleaned_sdrf_agrupado$scan_count <- sapply(strsplit(cleaned_sdrf_agrupado$Scans, ", "), length)
  
  num_scans <- max(cleaned_sdrf_agrupado$scan_count)
  samples <- sapply(1:num_scans, function(i) paste("Sample", i))
  
  output_path <- file.path(output_dir, "samples.txt")
  
  writeLines(samples, con = output_path)
  
  affy.set <- affy::ReadAffy(celfile.path = file.path(base_dir, "DatosExp", accesion_number))
                                                   
  eset <- affy::rma(affy.set)
  exprs_matrix <- as.data.frame(exprs(eset))
  
  exprs_matrix <- exprs_matrix[apply(exprs_matrix, 1, function(row) {
    all(!is.na(row) & trimws(row) != "")
  }), ]
  
  exprs_by_time <- list()
  cleaned_sdrf_agrupado <- as.data.frame(cleaned_sdrf_agrupado)
  colnames(cleaned_sdrf_agrupado)[1] <- "Tiempo"

  for (t in valores_tiempo) {
    scans_t <- cleaned_sdrf_agrupado %>%
      dplyr::filter(Tiempo == t) %>%
      dplyr::pull(Scans) %>%
      strsplit(", ") %>%
      unlist()
    
    exprs_by_time[[as.character(t)]] <- exprs_matrix %>%
      dplyr::select(all_of(scans_t))
  }
  
  resources <- list()
  for (t in names(exprs_by_time)) {
    write.table(exprs_by_time[[t]], 
                file = paste0(output_dir, "/exprs_time_", t, ".csv"),
                sep = ";",
                col.names = FALSE,
                row.names = FALSE,
                quote = FALSE)
    
    resources <- append(resources, paste0("exprs_time_", t, ".csv"))
  }
  
  
  ### ARCHIVO SÍMBOLOS DE GENES
  
  keys <- geneNames(affy.set)
  chip_name <- as.character(paste0(cleancdfname(affy.set@cdfName, addcdf = FALSE),".db"))
  
  if (!requireNamespace(chip_name, quietly = TRUE)) {
    BiocManager::install(chip_name, ask = FALSE)
  }
  
  library(chip_name, character.only = TRUE) 
  db <- get(chip_name)
  gene_symbols_mapIDs <- AnnotationDbi::mapIds(db, keys, column = c("SYMBOL"), 
                                               keytype = "PROBEID", multiVals = "first")
  
  output_path <- file.path(output_dir, "genes.txt")
  writeLines(gene_symbols_mapIDs, con = output_path)
  
  
  ##  MODIFICACION XML
  xml_path <- file.path(
    base_dir,
    "TrLab5", "resources", "resources.xml"
  )
  xml_file <- read_xml(xml_path)

  idf_path <- file.path(base_dir,"DatosExp", accesion_number, paste0(accesion_number, ".idf.txt"))
  idf_data <- read.delim(idf_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  idf_data <- t(idf_data)
  colnames(idf_data) <- as.character(idf_data[1, ])  
  idf_data <- as.data.frame(idf_data[-1, ])
  idf_data[idf_data == ""] <- NA
  
  datasets <- xml_find_all(xml_file, ".//dataset")
  num_datasets <- length(datasets)
  id <- num_datasets + 1
  organism <- grep("organism",colnames(sdrf_data), ignore.case = TRUE, value = TRUE)[1]
  
  exp_info <- c(idf_data$`Experiment Description`[1],
                length(keys),
                id,
                most_common_count,
                "150",
                length(cleaned_value_times),
                2,
                2,
                2,
                accesion_number,
                unique(sdrf_data[[organism]])[1],
                as.integer(num_scans),
                length(cleaned_value_times),
                "b"
  )
  
  if (tolower(exp_info[11]) == "mus musculus") {
    exp_info[11] <- "Mouse"
  }
  else if (tolower(exp_info[11]) == "homo sapiens"){
    exp_info[11] <- "Human"
  }
  
  dataset_name <- exp_info[10]
  
  existing <- xml_find_all(
    xml_file,
    sprintf(".//dataset[@name='%s']", dataset_name)
  )
  
  if (length(existing) == 0) {
    new_dataset <- xml_add_child(xml_file, "dataset", 
                                 description = exp_info[1],
                                 geneSize    = exp_info[2],
                                 id          = exp_info[10],
                                 maxC        = exp_info[4],
                                 maxG        = exp_info[5],
                                 maxT        = exp_info[6],
                                 minC        = exp_info[7],
                                 minG        = exp_info[8],
                                 minT        = exp_info[9],
                                 name        = dataset_name,
                                 organism    = exp_info[11],
                                 sampleSize  = exp_info[12],
                                 timeSize    = exp_info[13],
                                 type        = "b"
    )
    
    resources_node <- xml_add_child(new_dataset, "resources", separator = ";")
    for (resource in resources) {
      xml_add_child(resources_node, "resource", resource)
    }
    xml_add_child(new_dataset, "genes",   "genes.txt")
    xml_add_child(new_dataset, "samples", "samples.txt")
    xml_add_child(new_dataset, "times",   "time_points.txt")
    
    write_xml(xml_file, file.path(base_dir, "TrLab5/resources/resources.xml"))
  } else {
    message("El experimento '", dataset_name, "' ya está registrado en resources.XML")
  }
  
  ### GENERACIÓN ARCHIVO CONFIGURACIÓN .TRICFG
  out <- file.path(base_dir, "TrLab5/OutputTriGen", accesion_number)
  if (!dir.exists(out)) dir.create(out)
  
  output_trigen_dir <- paste0(out,"/",accesion_number)
  
  cfg_lines <- c(
    paste0("dataset = ", exp_info[10]),
    paste0("out = ", output_trigen_dir),
    paste0("fitness = msl"),
    paste0("N = 3"),
    paste0("G = 500"),
    paste0("I = 200"),
    
    
    paste0("Ale = 0.5"),
    paste0("Sel = 0.5"),
    paste0("Mut = 0.4"),
    paste0("Wf = 0.8"),
    paste0("Wg = 0.0"),
    paste0("Wc = 0.05"),
    paste0("Wt = 0.05"),
    paste0("WOg = 0.05"),
    paste0("WOc = 0.03"),
    paste0("Wot = 0.02"),
    
    paste0("minG = ", exp_info[8]),
    paste0("maxG = ", exp_info[5]),
    
    paste0("minC = ", exp_info[7]),
    paste0("maxC = ", exp_info[4]),
    
    paste0("minT = ", exp_info[9]),
    paste0("maxT = ", exp_info[6])
  )
  
  tricfg_path <- file.path(output_dir, "config_file.tricfg")
  writeLines(cfg_lines, tricfg_path)
  
  return(list(
    archivos = list.files(output_dir, full.names = TRUE, recursive = TRUE)
  ))
  
}


### FUNCIÓN PARA LANZAR EL ALGORITMO TRIGEN
runTriGen <- function(accesion_number, base_dir) {
  jar_path <- file.path(base_dir, "TrLab5/TriGenApp.jar")
  config_dir <- file.path(base_dir, "TrLab5/resources", accesion_number)
  config_path <- file.path(config_dir, "config_file.tricfg")
  
  print(base_dir)
  
  output_dir <- file.path(base_dir, "TrLab5/OutputTriGen", accesion_number)
  
  cmd <- paste("java -jar", shQuote(jar_path), shQuote(config_path))
  message("Ejecutando comando:\n", cmd)
  system(cmd)
  
  if (dir.exists(output_dir)) {
    archivos <- list.files(output_dir, full.names = TRUE, recursive = TRUE)
  } else {
    archivos <- character(0)
  }
  print(archivos)
  
  return(list(dir = output_dir, archivos = archivos))
}

''
### FUNCIÓN PARA POBLAR BUSCADOR CONJUNTOS DE ANOTACIÓN CONSUMIENDO API
getSupportedAnnotDatasets <- function() {
  resp <- httr::GET("http://pantherdb.org/services/oai/pantherdb/supportedannotdatasets")
  txt  <- httr::content(resp, as = "text", encoding = "UTF-8")
  dat  <- jsonlite::fromJSON(txt)
  
  annot_list <- dat$search$annotation_data_sets
  
  ids <- sapply(annot_list, function(x) {
    if (!is.null(x$annotation_data_type$id) && length(x$annotation_data_type$id) == 1) {
      x$annotation_data_type$id
    } else {
      NA_character_
    }
  }, USE.NAMES = FALSE)
  
  sort(unique(na.omit(ids)))
}
''


### FUNCIÓN CONSUMO API PARA ANÁLISIS SOBREREPRESENTACIÓN
runOverrepAnalysis <- function(genes, organism, annotDataSet, enrichmentTestType = "FISHER", correction = "FDR") {
  body <- list(
    geneInputList      = paste(genes, collapse = ","),
    organism           = organism,
    annotDataSet       = annotDataSet,
    enrichmentTestType = enrichmentTestType,
    correction         = correction
  )
  resp <- POST(
    url    = "https://pantherdb.org/services/oai/pantherdb/enrich/overrep",
    body   = body,
    encode = "form"
  )
  json <- httr::content(resp, as = "text", encoding = "UTF-8")
  dat  <- fromJSON(json)
 
  as.data.frame(dat$results$result)
}


### FUNCIÓN AUXILIAR OBTENER N TERMINOS MAS SOBREREPRESENTADOS
getTopTerms <- function(df, n = 10) {
  top_terms <- df %>%
    dplyr::arrange(desc(number_in_list)) %>%
    dplyr::slice_head(n = 10) %>%
    unnest_wider(term, names_sep = "_") %>%
    dplyr::rename(term_label = term_label, term_id = term_id)
}

### FUNCIÓN AUXILIAR PARA CREAR GRÁFICO
plotEnrichment <- function(top_terms, title = NULL) {
  p <- ggplot(top_terms, aes(x = "", y = number_in_list, fill = term_label,
                             text = paste0("Término: ", term_label, "<br>Genes: ", number_in_list))) +
    geom_bar(stat = "identity", width = 1) +
    labs(title = title, fill = "Término") +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  ggplotly(p, tooltip = "text")
}

### FUNCIÓN MOSTRAR METRICAS
showMetrics <- function(accesion_number, base_dir){
  metricas_path <- file.path(base_dir,
    "TrLab5/OutputTriGen",
    accesion_number,
    accesion_number,
    paste0(accesion_number,"_triq.csv")
  )
  
  csv <- read.csv(metricas_path, sep = ";")
  csv <- csv %>%
    dplyr::select(-BEST.SOLUTIONN, -BEST.TRIQN, -MEANN, -SDEVN)
}


### FUNCIÓN OBTENER GRÁFICAS VARIACIÓN NIVELES 
### DE EXPRESIÓN POR PUNTO DE TIEMPO SEGUN CONDICIÓN.
geneExpInteractive <- function(accesion_number, base_dir) {
  coords_dir <- file.path(base_dir,
    "TrLab5/OutputTriGen",
    accesion_number,
    accesion_number,
    "coordinates"
  )
  
  csvs <- list.files(coords_dir, pattern = "\\.csv$", full.names = TRUE)
  if (length(csvs) == 0) {
    warning("No hay archivos CSV en ", coords_dir)
    return(list())
  }
  
  all_plots <- list()
  
  for (path in csvs) {
    df <- read.csv(path, sep = ";", stringsAsFactors = FALSE)
    fname <- tools::file_path_sans_ext(basename(path))
    conds <- unique(df$s)
    
    plots_per_cond <- lapply(conds, function(cond) {
      dsub <- df[df$s == cond, ]
      plot_ly(dsub, x = ~t, y = ~el, color = ~factor(g),
              colors = "Set1", mode = "lines+markers", type = "scatter",
              hoverinfo = "text",
              text = ~paste("Gen:", g,
                            "<br>Tiempo:", t,
                            "<br>Expr:", el)
      ) %>%
        layout(
          title = paste0(fname, " — condición s=", cond),
          xaxis = list(title = "Tiempo (t)"),
          yaxis = list(title = "Expresión (el)"),
          legend = list(title = list(text = "<b>Gen</b>"))
        )
    })
    names(plots_per_cond) <- paste0(fname, "_s", conds)
    all_plots[[fname]] <- plots_per_cond
  }
  return(all_plots)
}


### FUNCIÓN CONSULTA GENES DE CADA TRICLUSTER SOLUCIÓN
queryGeneInfo <- function(accesion_number, organism, base_dir) {
  # Directorio con los .txt de genes
  geneLists_dir <- file.path(base_dir,
    "TrLab5/OutputTriGen",
    paste(accesion_number),
    accesion_number,
    "genes"
  )
  
  geneLists <- list.files(geneLists_dir, pattern = "\\.txt$", full.names = TRUE)
  
  results <- lapply(geneLists, function(fpath) {
    genes <- readLines(fpath, warn = FALSE)
    genes <- genes[genes != "" & !is.na(genes)]

    params <- list(
      organism      = organism,
      geneInputList = paste(genes, collapse = ",")
    )
    resp <- GET(
      url   = "https://pantherdb.org/services/oai/pantherdb/geneinfo",
      query = params,
      accept_json()
    )
    print(resp)
    
    if (resp$status_code != 200) {
      warning("HTTP ", resp$status_code, " en geneinfo para ", basename(fpath))
      return(NULL)
    }
    
    dat <- fromJSON(httr::content(resp, as = "text", encoding = "UTF-8"))
    
    if (!is.null(dat$search$mapped_genes$gene)) {
      df_genes <- as.data.frame(dat$search$mapped_genes$gene,
                                stringsAsFactors = FALSE)
    } else {
      df_genes <- data.frame()
    }
  })
  
  names(results) <- tools::file_path_sans_ext(basename(geneLists))
  return(results)
}


consultaExp <- function(term) {
  clean_term <- gsub(" ", "_", term)
  url <- paste0(
    "https://www.ebi.ac.uk/biostudies/api/v1/ArrayExpress/search?experimental_design=",
    URLencode(clean_term, reserved = TRUE),
    "&&link_type=GEO&&pageSize=100"
  )

  res <- httr::GET(url)
  httr::stop_for_status(res)
  raw <- httr::content(res, as="text", encoding="UTF-8")
  
  dat  <- fromJSON(raw, flatten = TRUE)
  df <- dat$hits[, c("accession", "title", "files")]

  return(df)
}








