########## Analisis diferencial de expresion RNAseq con EdgeR ############

# Salmon’s main output is its quantification file. This file is a plain-text, tab-separated file with a single header line (which names all of the columns). This file is named quant.sf and appears at the top-level of Salmon’s output directory. The columns appear in the following order:
#         Name, Length, EffectiveLength, TPM, NumReads
# Each subsequent row describes a single quantification record. The columns have the following interpretation.
  # Name — This is the name of the target transcript provided in the input transcript database (FASTA file).
  # Length — This is the length of the target transcript in nucleotides.
  # EffectiveLength — This is the computed effective length of the target transcript. It takes into account all factors being modeled that will effect the probability of sampling fragments from this transcript, including the fragment length distribution and sequence-specific and gc-fragment bias (if they are being modeled).
  # TPM — This is salmon’s estimate of the relative abundance of this transcript in units of Transcripts Per Million (TPM). TPM is the recommended relative abundance measure to use for downstream analysis.
  # NumReads — This is salmon’s estimate of the number of reads mapping to each transcript that was quantified. It is an “estimate” insofar as it is the expected number of reads that have originated from each transcript given the structure of the uniquely mapping and multi-mapping reads and the relative abundance estimates for each transcript.

# Lo primero que hay que hacer es instalar EdgeR si no lo tenemos instalado
BiocManager::install("edgeR")
BiocManager::install("tximport")

library(limma)
library(edgeR)
library(dplyr)
library(tidyr)
library(tximport)
library(readr)


# Cargamos el diseño experimental y la traducción de los transcritos a genes
sample_info <- read.csv("Design.csv") #sirve para construir rutas a los resultados
tx2gene <- read_tsv("Transcrito_a_Gen.tsv", col_names = FALSE) #lee el tsv indicando que el archivo no tiene cabecera
# este paso es esencial porque salmon quantifica a nivel de transcrio pero los analisis bio se hacen a nivel de gen, es el puente ente transcriptomica y genomica
# tx to gene (tx2gene) transcript to gene
colnames(tx2gene) <- c("TXNAME", "GENEID") # asigna nombres  por convencion, a las dos columnas de tx2gene con txname y geneid, sirve para que tximport pueda utilizarlo
# ya que necesita la primera columna de IDs de transcritos y la segunda de IDs de genes

# Definimos rutas a los ficheros de cuantificacion de Salmon
files <- file.path("quant_results", paste0(sample_info$Sample, "_quant"), "quant.sf") 
names(files) <- sample_info$Sample

# Leemos los datos de expresión con tximport, ojo, en design.csv hay que tener solo los que necesitamos o tximport fallara
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

# Exploración de la matriz
counts <- txi$counts
group <- factor(sample_info$Condition)
# Extrae la matriz de conteos por gen
# Esta matriz es:
#   La entrada directa a DESeq2, edgeR, etc.
#   La base de todo analisis diferencial


##### ----------------Edge R-----------------------

# ---------------------------------------
# Crear el objeto DGEList (guardar el conteo y asociarlo a las muestras)
# ---------------------------------------

y <- DGEList(counts=counts,group=group) # es el nucleo de edgeR aqui se guardan los conteos y se asocian las muestras con los grupos

# ---------------------------------------
# Filtrado de genes por expresion (esto se haria de normal pero nosotros no porque tenemos muy pocos (segun enunciado))
# ---------------------------------------
# keep <- filterByExpr(y) # guardamos los que mas expresion tienen
# y <- y[keep,,keep.lib.sizes=FALSE]

# ---------------------------------------
# Normalizacion de los datos
# ---------------------------------------
y <- calcNormFactors(y)
# utilizamos esta funcion porque lo que queremos comparar es la expresion de cada gen entre diferentes muestras
# tenemos que normalizar porque EdgeR no lo hace:
# corregimos los datos para que las diferencias entre muestras reflejen verdaderos camnios biologicos y no por problemas tecnicos:
# En RNAseq cada muestra puede variar porque:
    # 1. tamaño de libreria: algunas muestras tienen mas lecturas que otras
    # 2. sesgo de composicion, si unos genes estan super expresados, los genes bajos aprece que tengan demasiado poco
    # 3. otros factores: diferencuas de preparacion de libreria, eficiencia de secuenciacion, batch effects

# ---------------------------------------
# Definir el diseño estadistico: que vamos a comparar contra que
# ---------------------------------------
design <- model.matrix(~group) # esto modela expresion vs muestra, aqui se decide que comparacion biologica hacemos


# ---------------------------------------
# Estimar la dispersion
# ---------------------------------------
y <- estimateDisp(y, design)
  # edgeR modela los conteos con una binomial negativa y estima dispersion comun y dispersion por gen

# ---------------------------------------
# Ajustar el modelo y test estadistico
# ---------------------------------------
fit <- glmQLFit(y, design) 
qlf <- glmQLFTest(fit, coef = 2) # glmQFTestcompara los levels por la columna 2 (coef=2) es decir: Normopeso sera la base y vera si los genes se sobreexpresan o no en Sobrepeso

# Guardamos todos los genes
results <- topTags(qlf, n = Inf)$table
 # -> logFC es el cambio de expr -- Positivo → mas expresado en la condicion de interes. Negativo → menos expresado
 # -> pvalue es la significacion
 # -> FDR es lo que se reporta (false discovery rate) proporcion esperada de falsos positivos entre los genes significativos (Es un pval corregido)
    # FDR < 0.05 es el estandar cientifico

deg <- results[results$FDR < 0.05 & abs(results$logFC) > 1, ]
degFDR <- results[results$FDR < 0.05, ]


# Resumen EdgeR ----------
#  Conteos → DGEList
#  Filtrado → quitar ruido
#  Normalización → hacer muestras comparables
#  Modelo → biología experimental
#  Dispersión → variabilidad real
#  Test → genes diferenciales



#------------------------- Visualizacion de datos -----------------------------

# ---------------------------------------
# Definir umbrales de significancia
# ---------------------------------------
FDR_threshold <- 0.05
logFC_threshold <- 1 #(significa que la expresion se duplica o se reduce a la mitad)

status <- decideTests(qlf, p.value = FDR_threshold)


# ---------------------------------------
#Generar el grafico MA plot (el de edgeR)
# ---------------------------------------
# 1. Generar el gráfico
# values=c(1, -1) indica que queremos colorear los que suben (1) y bajan (-1)
# col=c("blue", "red") asigna azul a los que bajan y rojo a los que suben (puedes cambiarlo)
plotMD(qlf, status = status, values = c(1, -1), col = c("red", "blue"),
       main = "MA Plot: Normopeso vs Sobrepeso tipo 2",
       hl.cex = 0.8, legend=FALSE) # Tamaño de los puntos resaltados
# 2. Añadir leyenda (opcional pero recomendado)
legend("topleft", legend=c("Up-regulated", "Down-regulated", "Not Significant"),
       col=c("red", "blue", "black"), pch=16, cex=0.8, bty = "n")

# 3. (Opcional) Añadir líneas horizontales en logFC 0
abline(h=0, col="gray", lty=2)
abline(h=1, col="black", lty=2)
abline(h=-1, col="black", lty=2)

genes_candidatos <- results[abs(results$logFC) > 0.7, ]

ajuste_y <- ifelse(rownames(genes_candidatos) == "BDNF", 0.3, 0)
# 4. Pintamos las etiquetas limpias
text(x = genes_candidatos$logCPM, 
     y = genes_candidatos$logFC + ajuste_y, 
     labels = rownames(genes_candidatos), 
     cex = 0.7,   # Tamaño letra
     pos = 1,     # A la derecha del punto
     col = "black")





library(EnhancedVolcano)
library(ggplot2)
library(ggrepel)

# --- PASO 1: DEFINIR LOS COLORES (Crear el objeto 'keyvals') ---

# Primero ponemos todo en negro ('NS' = No Significativo)
keyvals <- rep('black', nrow(results))
names(keyvals) <- rep('NS', nrow(results))

# Pintamos de ROJO los que suben (logFC > 1 y FDR < 0.05)
# Asegúrate de que tus columnas se llaman logFC y FDR
keyvals[which(results$logFC > 1 & results$FDR < 0.05)] <- 'red'
names(keyvals)[which(results$logFC > 1 & results$FDR < 0.05)] <- 'Up-regulated'

# Pintamos de AZUL los que bajan (logFC < -1 y FDR < 0.05)
keyvals[which(results$logFC < -1 & results$FDR < 0.05)] <- 'blue'
names(keyvals)[which(results$logFC < -1 & results$FDR < 0.05)] <- 'Down-regulated'

EnhancedVolcano(results, lab = rownames(results),
                x = "logFC",
                y = "FDR",
                selectLab = rownames(results)[which(names(keyvals) %in% c('Up-regulated', 'Down-regulated'))], # Etiquetar solo los significativos
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~-Log[10]~ 'FDR'),
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 1.5,
                labSize = 4.5,
                colCustom = keyvals,
                colAlpha = 0.8,
                legendPosition = 'right',
                legendLabSize = 10,
                legendIconSize = 2,
                )



#-------PHEATMAP---------
library(pheatmap)
# 1. Calculamos los logCPM de TODAS tus muestras
# Usamos el objeto 'qlf' o 'y' (tu objeto DGEList original)
logCPM <- cpm(y, log=TRUE)

# 2. Seleccionamos solo los Mejores Genes
# No podemos pintar 20,000 genes, así que cogemos los Top 30 o 50 de tu tabla de resultados
# Ordenamos por FDR (los más significativos primero)
nombres_genes <- rownames(results)[order(results$FDR)]

# 3. Recortamos la matriz: Nos quedamos solo con las filas de esos genes top
matriz_heatmap <- logCPM[nombres_genes, ]

# 4. dibujamos el pheatmap
pheatmap(matriz_heatmap,
         scale = "row",           # ¡VITAL! Estandariza por filas para ver diferencias relativas
         cluster_rows = TRUE,     # Agrupa los genes que se comportan igual
         cluster_cols = TRUE,     # Agrupa los pacientes que se parecen entre sí
         cutree_rows = 4,  # Corta el árbol de genes en 2 bloques (Suben vs Bajan)
         cutree_cols = 3,
         clustering_method = "ward.D2",
         show_rownames = TRUE,    # Muestra nombres de genes
         show_colnames = TRUE,    # Muestra nombres de muestras
         main = "Genes Diferencialmente Expresados",
         fontsize_row = 8,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50) # Colores azul-blanco-rojo
)


