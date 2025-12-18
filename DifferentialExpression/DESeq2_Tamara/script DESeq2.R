
#### Script actividad grupal: Análisis de expresión diferencial con DESeq2 ####

install.packages("BiocManager")
BiocManager::install(c("tximport", "DESeq2", "readr"))

library(tximport)
library(DESeq2)
library(readr)


## Cargar el diseño experimental

design <- read_csv("data/Design.csv")

design


## Crear el vector con las rutas a quant.sf

files <- file.path(
  "data/quant_results",
  design$Sample,
  "quant.sf"
)

names(files) <- design$Sample

design$Condition <- gsub("/", "_", design$Condition)
design$Condition <- factor(design$Condition)


files


## Cargar la tabla Transcrito_a_Gen

tx2gene <- read_tsv(
  "data/Transcrito_a_Gen.tsv",
  col_names = c("TXNAME", "GENEID")
)

tx2gene$TXNAME <- sub("\\..*", "", tx2gene$TXNAME)

head(tx2gene)


## Importar las cuantificaciones con tximport

txi <- tximport(
  files,
  type = "salmon",
  tx2gene = tx2gene,
  ignoreTxVersion = TRUE
)


## Crear el objeto DESeq2

dds <- DESeqDataSetFromTximport(
  txi = txi,
  colData = design,
  design = ~ Condition
)

# Filtrado básico
dds <- dds[rowSums(counts(dds)) > 10, ]


## Ejecutar DESeq2
dds <- DESeq(dds)


## Comparación: Sobrepeso/Obeso2 vs Normopeso

res <- results(
  dds,
  contrast = c("Condition", "Sobrepeso_Obeso2", "Normopeso"),
  alpha = 0.05
)

res <- lfcShrink(
  dds,
  coef = "Condition_Sobrepeso_Obeso2_vs_Normopeso",
  type = "apeglm"
)


## Ver resultados 

summary(res)
head(res[order(res$padj), ])


## Genes significativos 
sig <- res[which(res$padj < 0.05 & abs(res$log2FoldChange) >= 1), ]
nrow(sig)



## Visualización


### Volcano Plot

# Preparar data frame para EnhancedVolcano
# Convertir a dataframe y añadir los nombres de los genes
res_df <- as.data.frame(res)
res_df$Gene <- rownames(res_df)

# Para volcano plot: columna con significancia ya creada
res_df$signif <- "No"
res_df$signif[res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1] <- "Sí"

head(res_df)

# Volcano plot con EnhancedVolcano

if(!require(EnhancedVolcano)){
  BiocManager::install("EnhancedVolcano")
  library(EnhancedVolcano)
}

EnhancedVolcano(
  res_df,
  lab = res_df$Gene,
  x = 'log2FoldChange',
  y = 'padj',        # padj para eje Y, más conservador
  pCutoff = 0.05,
  FCcutoff = 1,      # log2 fold change mínimo para marcar en volcano
  labSize = 3,
  axisLabSize = 10,
  pointSize = 2.5,
  colAlpha = 0.7,
  title = "Volcano plot: Obeso2 vs Normopeso",   # Título principal
  subtitle = "Genes diferencialmente expresados (p < 0.05)", # Subtítulo
  legendPosition = "bottom", # Posición de la leyenda
  titleLabSize = 14,
  subtitleLabSize = 12,
  legendLabSize = 10,
  legendIconSize = 3,
  drawConnectors = FALSE,   # Líneas conectando etiquetas
)

## Heatmap de los genes significativos

if(!require(pheatmap)){
  install.packages("pheatmap")
  library(pheatmap)
}

# Transformación VST (variance-stabilizing transformation)
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

# Extraer solo los genes significativos (padj < 0.05 y |LFC| >= 1)
sig_genes <- rownames(res)[which(res$padj < 0.05 & abs(res$log2FoldChange) >= 1)]

mat <- assay(vsd)[sig_genes, ]

# Escalar filas (genes)
mat_scaled <- t(scale(t(mat)))


# Nombres de las muestras
samples <- colnames(mat_scaled)

library(RColorBrewer)

# Tabla de anotaciones
annotation_col <- data.frame(
  Condition = design$Condition,
  Sexo = design$Sexo,
  Edad = design$Edad,
  row.names = design$Sample
)

# Colores correctos para pheatmap
annotation_colors <- list(
  Condition = c(
    Normopeso = "#A6CEE3",       
    Sobrepeso_Obeso2 = "#FB9A99"
  ),
  Sexo = c(
    Femenino = "#B2DF8A",
    Masculino = "#FDBF6F"
  ),
  Edad = colorRampPalette(c("#FFFFCC", "#FFCC66"))(length(unique(annotation_col$Edad))),
  Sample = setNames(
    colorRampPalette(brewer.pal(12, "Paired"))(nrow(annotation_col)),
    rownames(annotation_col)  # cada color asignado a cada muestra
  )
)


# Heatmap

pheatmap(
  mat_scaled,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 9,
  fontsize_col = 9,
  color = colorRampPalette(c("blue", "white", "red"))(50),  # expresión
  border_color = "grey80",
  main = "Heatmap de expresión de genes significativos"
)




## Boxplots de genes clave

# Top 5 genes más significativos
top_genes <- head(sig_genes[order(res$padj[sig_genes])], 5)

library(reshape2)

# Extraer datos de conteo normalizados (VST)
vst_mat <- assay(vsd)[top_genes, ]
vst_df <- melt(vst_mat)
colnames(vst_df) <- c("Gene", "Sample", "Expression")

# Añadir información de la condición
vst_df$Condition <- design$Condition[match(vst_df$Sample, design$Sample)]

# Boxplot con ggplot2
ggplot(vst_df, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_boxplot() +
  facet_wrap(~ Gene, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))



## MA Plot

plotMA(res, main = "MA Plot: Obeso2 vs Normopeso", alpha = 0.05)



