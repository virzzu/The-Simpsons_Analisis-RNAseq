
setwd("C:/Users/yanni/Desktop/MASTER/CUATRI 1/Secuenciación y ómicas/Actividad 2")

# Cargar librerías necesarias
BiocManager::install("EnhancedVolcano")
BiocManager::install("edger")

library(tximport) # Para importar datos de Salmon
library(edgeR) # Para análisis de expresión diferencial
library(readr) # Para lectura de archivos
library(ggplot2)  # Para gráficos
library(pheatmap) # Para heatmaps
library(RColorBrewer) # Para paletas de colores
library(EnhancedVolcano)# Para volcano plot mejorado
library(ggrepel)# Para etiquetas sin solapamiento
library(reshape2) # Para transformar datos


# PASO 1: CONFIGURACIÓN DE ARCHIVOS Y METADATOS

# Definir los archivos de cuantificación de Salmon (.sf)
samples <- c("bart_quant.sf", "lisa_quant.sf", "maggie_quant.sf", 
             "marge_quant.sf", "patty_quant.sf", "selma_quant.sf")

# Crear data frame con metadatos de las muestras
# Define qué muestra pertenece a qué grupo experimental
sample_names <- c("Bart", "Lisa", "Maggie", "Marge", "Patty", "Selma")
metadata <- data.frame(
  sample = sample_names,
  file = samples,
  condition = c("Normopeso", "Normopeso", "Normopeso",  # Bart, Lisa, Maggie
                "Obeso2", "Obeso2", "Obeso2"),           # Marge, Patty, Selma
  familia = c("Simpson", "Simpson", "Simpson",           
              "Simpson", "Bouvier", "Bouvier"),
  row.names = sample_names
)


# PASO 2: IMPORTAR DATOS DE SALMON CON TXIMPORT
txi <- tximport(samples, type="salmon", txOut=TRUE)

# Asignar nombres de muestras a las columnas de las matrices
colnames(txi$counts) <- sample_names      # Matriz de conteos
colnames(txi$abundance) <- sample_names   # Matriz de abundancias (TPM)
colnames(txi$length) <- sample_names      # Matriz de longitudes efectivas

# Extraer la matriz de conteos por gen
# Esta matriz es la entrada directa a edgeR y la base de todo análisis diferencial
counts <- txi$counts
group <- factor(metadata$condition, levels=c("Normopeso", "Obeso2"))



# PASO 3: CREAR OBJETO DGEList PARA edgeR

cts <- round(counts)
y <- DGEList(counts=cts, group=group)



# PASO 4: FILTRADO DE TRANSCRITOS CON BAJA EXPRESIÓN (NO!!!!!!!!!!!!!)

# Filtrar transcritos con expresión muy baja para quitar ruido. Mantener solo transcritos con >1 CPM (counts per million) en al menos 2 muestras. Esto reduce ruido y mejora el poder estadístico
keep <- rowSums(cpm(y) > 1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]


# PASO 5: NORMALIZACIÓN TMM

# Normalización por TMM (Trimmed Mean of M-values)
# EdgeR no normaliza automáticamente, debemos hacerlo explícitamente
# Corrige diferencias en:

y <- calcNormFactors(y, method="TMM")



# PASO 6: DISEÑO EXPERIMENTAL Y ESTIMACIÓN DE DISPERSIÓN

# Crear matriz de diseño para el modelo lineal
design <- model.matrix(~ group)
colnames(design) <- c("Intercept", "Obeso2_vs_Normopeso")

# Estimar dispersión (variabilidad) de los datos
y <- estimateDisp(y, design=design)


# PASO 7: ANÁLISIS DE EXPRESIÓN DIFERENCIAL

# Ajustar modelo con Quasi-Likelihood (más robusto que GLM estándar)
# glmQLFit es más conservador y controla mejor la tasa de falsos positivos
fit <- glmQLFit(y, design)

# Test de Quasi-Likelihood F-test
# coef=2 compara los levels por la columna 2: Obeso2 vs Normopeso
qlf <- glmQLFTest(fit, coef=2)

# PASO 8: OBTENER Y GUARDAR RESULTADOS

# Obtener TODOS los genes analizados
results <- topTags(qlf, n=Inf)$table

# Guardar todos los resultados
write.csv(results, "edgeR_results_all.csv", row.names=TRUE)

# Filtrar genes significativos con FDR < 0.05 y |logFC| > 1
deg <- results[results$FDR < 0.05 & abs(results$logFC) > 1, ]
write.csv(deg, "edgeR_results_significant_FC1.csv", row.names=TRUE)

# Filtrar genes significativos solo por FDR < 0.05 (sin umbral de fold change)
degFDR <- results[results$FDR < 0.05, ]
write.csv(degFDR, "edgeR_results_significant_FDR.csv", row.names=TRUE)

# Guardar metadatos
write.csv(metadata, "metadata.csv", row.names=TRUE)





# VISUALIZACIONES


# 1. MA PLOT
# Muestra expresión promedio (eje X) vs cambio de expresión (eje Y)
# Permite detectar sesgos relacionados con el nivel de expresión

# Definir umbrales de significancia
FDR_threshold <- 0.05
logFC_threshold <- 1  # Significa que la expresión se duplica o se reduce a la mitad

# Decidir qué genes son significativos
status <- decideTests(qlf, p.value=FDR_threshold)

# Crear MA plot con colores personalizados
# values=c(1, -1) indica que queremos colorear los que suben (1) y bajan (-1)
# col=c("red", "blue") asigna rojo a los que suben y azul a los que bajan
plotMD(qlf, status=status, values=c(1, -1), col=c("red", "blue"),
       main="MA Plot: Normopeso vs Obeso2",
       hl.cex=1.5, bg.cex=1, legend=FALSE)

# Añadir leyenda
legend("topleft", legend=c("Sobreexpresado", "Infraexpresado", "No significativo"),
       col=c("red", "blue", "black"), pch=16, cex=0.8, bty="n")

# Añadir líneas de referencia
abline(h=0, col="gray", lty=2) # Línea de sin cambio
abline(h=1, col="black", lty=2) # Umbral superior
abline(h=-1, col="black", lty=2)  # Umbral inferior

# Etiquetar genes candidatos (con logFC > 0.7)
genes_candidatos <- results[abs(results$logFC) > 0.7, ]
if(nrow(genes_candidatos) > 0) {
  text(x=genes_candidatos$logCPM, 
       y=genes_candidatos$logFC, 
       labels=rownames(genes_candidatos), 
       cex=0.8, pos=1, col="black")
}

# 2. VOLCANO PLOT 
# Muestra la relación entre cambio de expresión (logFC) y significancia (FDR)

# Crear vector de colores personalizados
# Primero ponemos todo en gris (No Significativo)
keyvals <- rep('darkgrey', nrow(results))
names(keyvals) <- rep('No significativo', nrow(results))

# Rojo los que aumentan (logFC > 1 y FDR < 0.05)
keyvals[which(results$logFC > 1 & results$FDR < 0.05)] <- 'red'
names(keyvals)[which(results$logFC > 1 & results$FDR < 0.05)] <- 'Sobreexpresados'

# Azul los que bajan (logFC < -1 y FDR < 0.05)
keyvals[which(results$logFC < -1 & results$FDR < 0.05)] <- 'blue'
names(keyvals)[which(results$logFC < -1 & results$FDR < 0.05)] <- 'Infraexpresados'

# Identificar genes para etiquetar
genes_significativos <- rownames(results)[which(names(keyvals) %in% c('Sobreexpresados', 'Infraexpresados'))]
genes_limite <- rownames(results)[which(results$FDR < 0.05)]

# Crear volcano plot 
EnhancedVolcano(results,
                title="Expresión diferencial génica:\nNormopeso vs Obeso2",
                lab=rownames(results),
                x="logFC",
                y="FDR",
                ylim=c(-0.5, max(-log10(results$FDR), na.rm=TRUE) + 1),
                xlim=c(min(results$logFC, na.rm=TRUE) - 0.5, max(results$logFC, na.rm=TRUE) + 0.5),
                selectLab=unique(c(genes_significativos, genes_limite)),
                xlab=bquote(~Log[2]~ 'fold change'),
                ylab=bquote(~-Log[10]~ "FDR"),
                pCutoff=0.05,
                FCcutoff=1,
                pointSize=2.5,
                labSize=5.5,
                colCustom=keyvals,
                colAlpha=0.8,
                legendPosition='right',
                legendLabSize=10,
                legendIconSize=2,
                drawConnectors=TRUE,
                widthConnectors=0.01,
                colConnectors="black",
                arrowheads=FALSE)

# 3. HEATMAP
# Muestra patrones de expresión de los genes más diferencialmente expresados
# Agrupa muestras y genes según similitud

# Calcular logCPM de todas las muestras
logCPM <- cpm(y, log=TRUE)

# Seleccionar los top genes (ordenados por FDR, los más significativos primero)
# No podemos pintar 20,000 genes, así que cogemos los Top 30-50
n_top <- min(50, nrow(results))
nombres_genes <- rownames(results)[order(results$FDR)][1:n_top]

# Recortar la matriz: solo las filas de esos genes top
matriz_heatmap <- logCPM[nombres_genes, ]

# Crear tabla de anotaciones (leyenda de las columnas)
anotaciones_muestras <- data.frame(Grupo=y$samples$group)
rownames(anotaciones_muestras) <- colnames(matriz_heatmap)

# Definir colores para los grupos
colores_anotacion <- list(
  Grupo=c(Normopeso="pink", Obeso2="skyblue")
)

# Crear y guardar heatmap
pheatmap(matriz_heatmap,
         scale="row",  # Normalizar por filas (cada gen tiene media 0 y sd 1)
         annotation_col=anotaciones_muestras,  # Añadir barra de grupos
         annotation_colors=colores_anotacion,  # Colores de la barra
         cluster_rows=TRUE,
         cluster_cols=TRUE,
         cutree_rows=4,  # Cortar el árbol de genes en 4 bloques
         cutree_cols=2,  # Cortar el árbol de muestras en 2 bloques
         clustering_method="ward.D2",
         main="Expresión génica: Normopeso vs Obeso2",
         color=colorRampPalette(c("blue", "white", "red"))(100),
         filename="heatmap.png",
         width=10,
         height=12)
