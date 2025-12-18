################################### PRUEBA DEFINITI??? ##########################
setwd("~/Ómicas/Activ grupal")
library(DESeq2)
library(tximport)
library(EnhancedVolcano)
library(pheatmap)
library(ggplot2)
library(reshape2)


files <- c(
  BartSimpson = "C:/Users/vdlhh/OneDrive/Documentos/Ómicas/Activ grupal/BartSimpson_quant.sf",
  LisaSimpson = "C:/Users/vdlhh/OneDrive/Documentos/Ómicas/Activ grupal/LisaSimpson_quant.sf",
  MaggieSimpson = "C:/Users/vdlhh/OneDrive/Documentos/Ómicas/Activ grupal/MaggieSimpson_quant.sf",
  MargeSimpson = "C:/Users/vdlhh/OneDrive/Documentos/Ómicas/Activ grupal/MargeSimpson_quant.sf",
  PattyBouvier = "C:/Users/vdlhh/OneDrive/Documentos/Ómicas/Activ grupal/PattyBouvier_quant.sf",
  SelmaBouvier = "C:/Users/vdlhh/OneDrive/Documentos/Ómicas/Activ grupal/SelmaBouvier_quant.sf"
)

tx2gene <- read.delim(
  "Transcrito_a_Gen.tsv",
  header = TRUE,
  stringsAsFactors = FALSE
)

txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
counts <- txi$counts
head(counts)

coldata <- read.csv(
  "metadatos_actividad_grupal.csv",
  row.names = 1
)

coldata <- data.frame(
  row.names = colnames(counts),
  condition = c("normopeso", "normopeso", "normopeso", "obeso2", "obeso2", "obeso2")
)

stopifnot(all(rownames(coldata) == colnames(counts)))

tx2gene <- read.table(
  "Transcrito_a_Gen.tsv",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)


dds <- DESeqDataSetFromTximport(
  txi = txi,
  colData = coldata,
  design = ~ condition
)

dds <- DESeq(dds)
res <- results(dds)
head(res_df)
head(res)

plotMA(res)

res_df <- as.data.frame(res)
res_df <- res_df[!is.na(res_df$padj), ]
res_df$gene <- rownames(res_df)
res_df$signif <- res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1

library(ggplot2)

ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = signif), alpha = 0.8, size = 2) +   # TODOS los genes
  geom_text(
    data = subset(res_df, signif), 
    aes(label = gene),
    size = 3,
    vjust = -0.5,
    check_overlap = TRUE
  ) +
  coord_cartesian(xlim = c(-5, 5)) +
  scale_color_manual(
    values = c("grey40", "firebrick3")   
  ) +
  theme_minimal() +
  labs(
    x = "log2 Fold Change",
    y = "-log10(pvalue)",
    color = "Significativo"
  )

volcano_plot(as.data.frame(res))

if(!require(pheatmap)){
  install.packages("pheatmap")
  library(pheatmap)
}
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
mat <- assay(vsd)
mat <- mat - rowMeans(mat)

pheatmap(mat,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = as.data.frame(colData(dds)),
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 3,
         fontsize_col = 8,
         color = colorRampPalette(c("blue", "white", "red"))(50))

vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
mat <- assay(vsd)
plot_df <- melt(mat)
colnames(plot_df) <- c("gene", "sample", "expression")
plot_df$condition <- colData(dds)$condition[
  match(plot_df$sample, rownames(colData(dds)))
]

ggplot(plot_df, aes(x = sample, y = expression, fill = condition)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
  theme_minimal() +
  labs(
    x = "Muestra",
    y = "Expresión (VST)",
    fill = "Condición"
  )