#Autor: Samuel Pintos


setwd("~/BIOINFORMÁTICA_R/Secuenciación/ACT2/Quants")
# 1.1. Cargar paquetes necesarios

# 1.1. Instalar BiocManager (Si es necesario)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# 1.2. Instalar paquetes de análisis y anotación
# ¡IMPORTANTE! Los nombres de los paquetes van ENTRE COMILLAS
BiocManager::install(c("tximport", "rtracklayer", "limma", "edgeR", "AnnotationDbi", "org.Hs.eg.db"), force = TRUE)

# 1.3. CARGAR PAQUETES (¡Ejecutar CADA VEZ que inicies R!)
library(tximport)
library(rtracklayer) 
library(limma)
library(edgeR) 
library(AnnotationDbi) # Ahora sí funcionará
library(org.Hs.eg.db) # Anotación para Homo sapiens


# ASUME que todos tus archivos .sf están en la misma carpeta o que la ruta es correcta.
# Por simplicidad, usaremos el directorio de trabajo de R (wd) o una ruta directa.
# Reemplaza "ruta/a/tus/quants/" con la ruta real si no usas el wd.
# Nombres de tus muestras
sample_names <- c("Bart", "Patty", "Selma", "Lisa", "Marge", "Maggie")
# Nombres de los archivos completos (asumiendo que están en el directorio de trabajo o una ruta única)
files <- paste0(sample_names, "_quant.sf")
names(files) <- sample_names
# 1.3. Definir Metadatos (Diseño Experimental)
# **ESTO ES UNA ASUNCIÓN. DEBES CAMBIAR LOS GRUPOS SEGÚN TU EXPERIMENTO REAL**
# Ejemplo: Bart, Patty, Selma son "Control"; Lisa, Marge, Maggie son "Tratado".
targets <- data.frame(
  SampleID = sample_names,
  File = files,
  Condition = factor(c("Normopeso", "Obeso2", "Obeso2",
                       "Normopeso", "Obeso2", "Normopeso")),
  row.names = sample_names
)
print(targets)

cat("\nGenerando mapeo Transcripto a Gen (tx2gene)... Esto puede tomar unos segundos.\n")

# 3.1. Obtener todos los IDs de transcripto del primer archivo (ejemplo)
first_quant <- read.delim(files[1], header = TRUE, stringsAsFactors = FALSE)
tx_ids_with_version <- first_quant$Name

# 3.2. Quitar la versión (ej: .1, .2) del ID para mapear con la base de datos
tx_ids_no_version <- gsub("\\..*","", tx_ids_with_version)

# 3.3. Mapear RefSeq Accessions (REFSEQ) a Entrez Gene IDs (ENTREZID)
# Utilizamos org.Hs.eg.db para realizar el mapeo de IDs
tx_to_eg <- AnnotationDbi::mapIds(
  x = org.Hs.eg.db,
  keys = tx_ids_no_version,
  column = "ENTREZID", # Columna destino (Entrez Gene ID)
  keytype = "REFSEQ",  # Tipo de IDs de entrada (RefSeq Accession)
  multiVals = "first"  # Manejo de IDs duplicados
)

# 3.4. Crear la tabla 'tx2gene' requerida por tximport
tx2gene <- data.frame(
  TXNAME = tx_ids_with_version, 
  GENEID = tx_to_eg, 
  stringsAsFactors = FALSE
)
# Limpiar entradas sin mapeo (NA) y duplicados
tx2gene <- na.omit(tx2gene)
tx2gene <- unique(tx2gene)

# ----------------------------------------------------------------------
# 4. IMPORTAR Y AGREGAR CUENTAS con TXIMPORT
# ----------------------------------------------------------------------

cat("Importando cuantificaciones y agregando a nivel de gen...\n")
txi <- tximport(files, 
                type = "salmon", 
                tx2gene = tx2gene,
                txOut = FALSE, 
                countsFromAbundance = "lengthScaledTPM") 

# Crear la matriz de cuentas requerida para DGEList
counts_matrix <- round(txi$counts)
cat("\n¡Objeto 'counts_matrix' creado! Listo para el análisis.\n")

# Crear el objeto DGEList (Lista de Expresión Diferencial de Genes)
# El objeto guarda las cuentas, las librerías totales (column sums) y la información de grupo.
d0 <- DGEList(counts=counts_matrix, group=targets$Condition, lib.size = colSums(counts_matrix)) 

View(d0)
# ------------------------------------
# 3. Filtrado de baja expresión (Mejora la potencia estadística)
# ------------------------------------

# Se filtran los genes que son muy poco expresados (ruido)
# Mantenemos genes que tienen un CPM (counts-per-million) > 0.5 en al menos 3 muestras 
# (el tamaño del grupo más pequeño que definimos).
keep <- filterByExpr(d0, min.count = 0.5, min.total.count = 10, large.n = 3)
d <- d0[keep,,keep.lib.sizes=FALSE]

cat(paste0("Número de genes antes del filtrado: ", nrow(d0$counts), "\n"))
cat(paste0("Número de genes después del filtrado: ", nrow(d$counts), "\n"))


# ------------------------------------
# 4. Normalización de Librería
# ------------------------------------

# Se calculan los factores de normalización (TMM: Trimmed Mean of M-values)
d <- calcNormFactors(d, method="TMM")

# ------------------------------------
# 5. Definición del Modelo y Voom
# ------------------------------------

# 5.1 Crear la Matriz de Diseño (similar a la matriz de diseño en DESeq2)
# Esto define la relación entre los grupos que queremos comparar.
design <- model.matrix(~0 + targets$Condition)
colnames(design) <- levels(targets$Condition) # Nombrar las columnas con los grupos
rownames(design) <- rownames(targets)
print(design)

# 5.2 Aplicar la Transformación voom
# Transforma las cuentas, estima la relación media-varianza y asigna pesos.
v <- voom(d, design, plot=TRUE) 
# 

# ------------------------------------
# 6. Ajuste del Modelo Lineal y Resultados
# ------------------------------------

# 6.1 Ajustar el modelo lineal (Fit)
fit <- lmFit(v, design)


# 1. Crear la matriz de diseño (sin cambiar los nombres)
design_original <- model.matrix(~0 + targets$Condition) 

# 2. Renombrar la matriz de diseño antes de usarla en voom
# ¡Asegúrese de que los nombres de sus niveles sean correctos!
colnames(design_original) <- levels(targets$Condition) 
rownames(design_original) <- rownames(targets)
design <- design_original # Usamos la matriz renombrada

# ... (Continuar con voom(d, design, plot=TRUE), lmFit(v, design), etc.)

# 3. Definir y Contrastar Comparaciones utilizando los nombres exactos:
contrast.matrix <- makeContrasts(Normopeso_vs_Obeso2 = Normopeso - Obeso2, 
                                 levels=design) # Los niveles deben coincidir con colnames(design)


# 6.3 Aplicar la Matriz de Contrastes
fit2 <- contrasts.fit(fit, contrast.matrix)

# 6.4 Aplicar ajuste Bayesiano (Empirical Bayes, el método central de limma)
fit2 <- eBayes(fit2)


# ------------------------------------
# 7. Obtener y Guardar Resultados
# ------------------------------------

# 7.1. Obtener la tabla de resultados para la comparación definida
# 'number=Inf' asegura que obtienes todos los genes, no solo los primeros 10.
results <- topTable(fit2, coef="Normopeso_vs_Obeso2", number=Inf, sort.by="p")

cat("\n--- Resultados del Análisis (Top 10 Genes más significativos) ---\n")
print(head(results, 10))

# 7.2. Identificar los genes diferencialmente expresados (DEG)
# Criterio común:
# - p-valor ajustado (adj.P.Val o FDR) < 0.05
# - Cambio de pliegue (logFC) mayor a 1 (o menor a -1)
deg_genes <- subset(results, adj.P.Val < 0.05 & abs(logFC) > 1)

cat(paste0("\nNúmero de genes diferencialmente expresados (FDR < 0.05, |logFC| > 1): ", nrow(deg_genes), "\n"))


# 7.3. Guardar resultados en formato CSV
write.csv(results, "limma_voom_Resultados_Completos.csv", row.names = TRUE)
write.csv(deg_genes, "limma_voom_Genes_Diferenciales.csv", row.names = TRUE)

cat("\nAnálisis finalizado con éxito. \n")
cat("Los resultados completos están en 'limma_voom_Resultados_Completos.csv'.\n")
cat("Los genes significativos (DEG) están en 'limma_voom_Genes_Diferenciales.csv'.\n")





# ----------------------------------------------------------------------
# 3. IDENTIFICAR NOMBRES DE GENES
# ----------------------------------------------------------------------

library(AnnotationDbi)
library(org.Hs.eg.db)

# 1. Creamos un vector con los Entrez IDs que queremos buscar
# (Usamos los nombres de fila de la tabla de resultados)
gene_ids_to_search <- rownames(results)[1:5] 

# 2. Mapeamos los IDs a Símbolos de Gen (Gene Symbol) y Nombre Completo
gene_symbols <- mapIds(org.Hs.eg.db, 
                       keys=gene_ids_to_search, 
                       column="SYMBOL", 
                       keytype="ENTREZID",
                       multiVals="first")

gene_names <- mapIds(org.Hs.eg.db, 
                     keys=gene_ids_to_search, 
                     column="GENENAME", 
                     keytype="ENTREZID",
                     multiVals="first")

# 3. Añadimos esta información a la tabla de resultados
top5_results <- head(results, 5)

top5_results$Symbol <- gene_symbols[rownames(top5_results)]
top5_results$GeneName <- gene_names[rownames(top5_results)]

# Reordenamos las columnas para que los nombres aparezcan primero
top5_final <- top5_results[, c("Symbol", "GeneName", "logFC", "AveExpr", "adj.P.Val")]

cat("\n--- Top 5 Genes con Nombres y Dirección del Cambio ---\n")
print(top5_final)

# Si quieres buscar un solo ID, por ejemplo el 23304:
# select(org.Hs.eg.db, keys="23304", columns=c("SYMBOL", "GENENAME"), keytype="ENTREZID")
