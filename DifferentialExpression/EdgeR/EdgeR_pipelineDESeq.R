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
# Extrae la matriz de conteos por gen:
# Esta matriz es:
#   La entrada directa a DESeq2, edgeR, etc.
#   La base de todo analisis diferencial


##### ----------------Edge R-----------------------
y <- DGEList(counts=counts,group=group) # es el nucleo de edgeR aqui se guardan los conteos y se asocian las muestras con los grupos

# Filtrado de genes por expresion
keep <- filterByExpr(y) # guardamos los que mas expresion tiene
y <- y[keep,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
