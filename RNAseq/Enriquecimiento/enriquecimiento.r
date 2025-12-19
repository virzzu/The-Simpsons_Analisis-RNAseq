#--------------------------------- Enriquecimiento ------------------------------------
# Autora: Virginia Garcia-Loygorri

#########################       Enriquecimiento       ##############################
# Carga de datos de edgeR para tratar con ellos
genes_DEG <- readRDS("../Expresion_diferencial/EdgeR/degFDR.rds")

# Cargamos librerias para enriquecimiento ! Con acceso a internet
library(gprofiler2)
library(plotly)


gost_data <- gost(query = list(" " = row.names(genes_DEG[genes_DEG$FDR < 0.05 & abs(genes_DEG$logFC) >= 1, ])), 
                 organism = "hsapiens", 
                 ordered_query = FALSE,
                 correction_method = "fdr",
                 sources = c("GO:BP", "KEGG", "REAC", "WP", "HPA", "HP"),
                 significant = TRUE)

# Generar el plot interactivo
p_interactivo <- gostplot(gost_data, interactive = TRUE)

# Personalizacion del plot
layout(p_interactivo,
       title = "\n\nAnálisis de enriquecimiento funcional de los genes expresados diferencialmente",
       xaxis = list(title = "Bases de Datos"))


#########################       Tabla de enrich       ##############################

gost_table <- gost_data$result[, c("source", "term_id", "term_name", "p_value", "intersection_size")]
#Cogemos los valores del enriquecimiento de gost que nos interesan

library(reactable) #libreria para tabla interactiva

reactable(
  gost_table,
  columns = list(
    source = colDef(name = "Base de Datos", width = 100, align = "center"),
    term_id = colDef(name = "ID", width = 120, align = "center"),
    term_name = colDef(name = "Descripción del Término", minWidth = 200, align = "center"),
    p_value = colDef(name = "P-Value Adj", align = "center", cell = function(value) {
      formatC(value, format = "e", digits = 3) }
      ),
    intersection_size = colDef(name = "Nº Genes", align = "center")
    ),
  pagination = TRUE,
  showPagination = TRUE,
  showPageInfo = TRUE,
  sortable = TRUE,
  searchable = TRUE,      # Barra de búsqueda
  striped = TRUE,         # Filas a rayas
  highlight = TRUE,       # Resalta al pasar el ratón
  bordered = TRUE,
  defaultPageSize = 15,   # Paginación
  theme = reactableTheme(
    style = list(fontFamily = "Helvetica, sans_serif",
    fontSize = "14px"),
    borderColor = "black",
    stripedColor = "#D5DBDE",
    highlightColor = "#C4E4F5",
    cellPadding = "8px 12px",
    headerStyle = list(
      fontWeight = "bold",
      backgroundColor = "#f0f0f0"
    )
  )
)
