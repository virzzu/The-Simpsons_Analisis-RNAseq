# ==============================================================================
# SCRIPT FINAL: ANÁLISIS DE EXPRESIÓN DIFERENCIAL (RNA-SEQ)
# Visualización en RStudio + Exportación automática a PNG
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. CARGA DE LIBRERÍAS
# ------------------------------------------------------------------------------
rm(list=ls()) 
library(edgeR)
library(limma)
library(ggplot2)
library(EnhancedVolcano)
library(pheatmap)
library(tximport)
library(grid)

# Establecer directorio de trabajo
setwd("C:/Users/ritap/OneDrive/Documents/Master Bioinformatica - UNIR/Sequenciació i Òmiques de 1a Generació")

# ------------------------------------------------------------------------------
# 2. CARGA DE DATOS Y METADATOS
# ------------------------------------------------------------------------------
metadata <- read.csv("metadatos_actividad_grupal(in).csv", row.names = 1)
metadata$Condition <- factor(metadata$Condition, levels = c("Normopeso", "Obeso2"))

ruta_principal <- getwd()
carpetas_quant <- list.dirs(ruta_principal, full.names = TRUE, recursive = FALSE)
carpetas_quant <- carpetas_quant[grep("_quant$", carpetas_quant)]

files <- file.path(carpetas_quant, "quant.sf")
names(files) <- gsub("_quant$", "", basename(carpetas_quant))
files <- files[rownames(metadata)]

# ------------------------------------------------------------------------------
# 3. IMPORTACIÓN (Salmon -> Gene Level)
# ------------------------------------------------------------------------------
tx2gene <- read.delim("Transcrito_a_Gen.tsv", header = FALSE)
colnames(tx2gene) <- c("TXNAME", "GENEID")
tx2gene$TXNAME <- gsub("\\..*$", "", tx2gene$TXNAME)

txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)

# ------------------------------------------------------------------------------
# 4. PROCESAMIENTO LIMMA-VOOM
# ------------------------------------------------------------------------------
dge <- DGEList(txi$counts)
dge <- calcNormFactors(dge)
keep <- filterByExpr(dge, group = metadata$Condition)
dge <- dge[keep, , keep.lib.sizes = FALSE]

design <- model.matrix(~0 + Condition, data = metadata)
colnames(design) <- levels(metadata$Condition)

# --- PLOT 1: VOOM ---
# 1. Visualizar en RStudio
v <- voom(dge, design, plot = TRUE) 

# 2. Guardar en PNG
png("1_Voom_Plot.png", width = 1000, height = 800, res = 150)
voom(dge, design, plot = TRUE)
dev.off()

# ------------------------------------------------------------------------------
# 5. ANÁLISIS ESTADÍSTICO
# ------------------------------------------------------------------------------
fit <- lmFit(v, design)
cont_matrix <- makeContrasts(Obeso_vs_Normal = Obeso2 - Normopeso, levels = design)
fit2 <- contrasts.fit(fit, cont_matrix)
fit2 <- eBayes(fit2)

res_limma <- topTable(fit2, coef = "Obeso_vs_Normal", number = Inf, sort.by = "P")
write.csv(res_limma, "Resultados_Limma_Obeso_vs_Normal.csv")

# ------------------------------------------------------------------------------
# 6. GENERACIÓN Y EXPORTACIÓN DE GRÁFICOS FINALES
# ------------------------------------------------------------------------------

# --- PLOT 2: VOLCANO ---
p_volcano <- EnhancedVolcano(res_limma,
                             lab = rownames(res_limma),
                             x = 'logFC', y = 'adj.P.Val',
                             pCutoff = 0.05, FCcutoff = 1.0,
                             title = 'Obeso2 vs Normopeso',
                             subtitle = 'FDR < 0.05 y |logFC| > 1',
                             legendPosition = 'bottom')

# 1. Mostrar en el panel Plots
print(p_volcano) 

# 2. Guardar en PNG
ggsave("2_Volcano_Plot.png", plot = p_volcano, width = 8, height = 7, dpi = 300)

# --- PLOT 3: HEATMAP ---
set.seed(1995)
top_genes <- rownames(res_limma)[1:30]
counts_log <- v$E[top_genes, ]
colnames(counts_log) <- rownames(metadata)

df_annot <- data.frame(Condicion = metadata$Condition)
rownames(df_annot) <- rownames(metadata)

# 1. Mostrar en RStudio (Viewer)
pheatmap(counts_log, 
         annotation_col = df_annot, 
         scale = "row", 
         main = "Top 30 Genes Diferenciales",
         cluster_cols = TRUE,
         cutree_rows = 3, cutree_cols = 2,
         color = colorRampPalette(c("navyblue", "white", "firebrick3"))(100))

# 2. Guardar en PNG (usando el argumento filename)
pheatmap(counts_log, 
         annotation_col = df_annot, 
         scale = "row", 
         main = "Top 30 Genes Diferenciales",
         cluster_cols = TRUE,
         cutree_rows = 3, cutree_cols = 2,
         color = colorRampPalette(c("navyblue", "white", "firebrick3"))(100),
         filename = "3_Heatmap.png", 
         width = 8, height = 10)

cat("\n>>> ¡PROCESO COMPLETADO! <<<\n")
cat("1. Revisa los gráficos en los paneles de RStudio.\n")
cat("2. En tu carpeta de trabajo tienes: 1_Voom_Plot.png, 2_Volcano_Plot.png, 3_Heatmap.png y el CSV de resultados.")