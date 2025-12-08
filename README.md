# The Simpsons Analisis-RNAseq
## Análisis de Expresión Diferencial de Genes Relacionados con la Obesidad (RNA-seq) en personajes de Los Simpson

**Asignatura:** Secuenciación y Ómicas de próxima generación
**Estudiantes:** 

## 1\. Descripción del Proyecto

Este proyecto consiste en un flujo de trabajo bioinformático para el análisis de datos de RNA-seq simulados. El objetivo es identificar genes diferencialmente expresados relacionados con la obesidad comparando perfiles genéticos de personajes del universo de *Los Simpson*.

El estudio busca encontrar marcadores moleculares y rutas metabólicas alteradas que expliquen el fenotipo de obesidad en el grupo asignado.

## 2\. Diseño Experimental

Este análisis se centra exclusivamente en la siguiente comparativa:

  * **Grupo Caso (Obesidad Tipo 2):** Fenotipo con sobrepeso/obesidad (Familia Bouvier/Simpson).
      * *Muestras:* Marge Simpson, Patty Bouvier, Selma Bouvier.
  * **Grupo Control (Normopeso):** Fenotipo de peso normal.
      * *Muestras:* Bart Simpson, Lisa Simpson, Maggie Simpson.

> **Nota:** Las muestras correspondientes al "Grupo Obeso 1" (Homer y Abraham) han sido excluidas de este análisis específico según la asignación del grupo de trabajo.

## 3\. Estructura de Datos

El proyecto asume la siguiente organización de directorios:

```text
.
├── data/
│   ├── raw/            # Archivos FASTQ originales (Marge, Patty, Selma, Bart, Lisa, Maggie)
│   ├── trimmed/        # Archivos FASTQ procesados tras el control de calidad
│   └── reference/      # Genoma de referencia y anotaciones (.fa, .gtf)
├── results/
│   ├── qc/             # Informes de calidad (FastQC, MultiQC)
│   ├── alignment/      # Archivos de alineamiento (.bam)
│   ├── counts/         # Matriz de conteos (counts.txt)
│   └── figures/        # Gráficos generados (Volcano plots, Heatmaps)
├── scripts/            # Scripts de Bash o R utilizados en el pipeline
└── README.md           # Documentación del proyecto
```

## 4\. Flujo de Trabajo (Pipeline)

### Paso 1: Control de Calidad (QC)

Evaluación de la calidad de las lecturas crudas para detectar adaptadores, bases de baja calidad o contaminación.

  * **Herramienta:** FastQC / MultiQC.
  * **Criterio:** Se revisaron métricas como *Per base sequence quality* y *Adapter content*.

### Paso 2: Limpieza y Recorte (Trimming)

Eliminación de adaptadores y recorte de bases de baja calidad en los extremos de las lecturas.

  * **Herramienta:** Fastp.
  * **Parámetros:** Eliminación de lecturas con longitud \< 35 pb y calidad media \< 30.

### Paso 3: Alineamiento (Mapping)

Mapeo de las lecturas procesadas contra el genoma de referencia simulado.

  * **Herramienta:** HISAT2 / STAR / Bowtie2.
  * **Output:** Archivos BAM ordenados y con índice.

### Paso 4: Cuantificación

Conteo de lecturas mapeadas en cada gen (feature).

  * **Herramienta:** featureCounts (Subread package) o HTSeq.
  * **Resultado:** Matriz de conteos crudos (Genes en filas, Muestras en columnas).

### Paso 5: Análisis de Expresión Diferencial

Análisis estadístico en R para identificar cambios significativos entre el **Grupo Obeso 2** y el **Grupo Normopeso**.

  * **Herramientas:** Paquete R (DESeq2 o edgeR).
  * **Diseño:** `~ condición` (donde condición es Obeso2 vs Normopeso).
  * **Criterios de significancia:**
      * |Log2 Fold Change| \> 1 (o el umbral que decidas, ej. 0.58 para 1.5 veces de cambio).
      * P-value ajustado (FDR) \< 0.05.

## 5\. Resultados Clave

Los resultados principales se encuentran en la carpeta `results/figures` e incluyen:

1.  **Volcano Plot:** Visualización global de la distribución de genes (significancia vs. magnitud de cambio).
2.  **Heatmap:** Mapa de calor de los 50 genes más variables, mostrando el agrupamiento (clustering) entre las hermanas Bouvier/Marge y los niños Simpson.
3.  **Tabla de DEGs:** Listado de genes candidatos con su anotación funcional.

## 6\. Requisitos del Sistema

Para reproducir este análisis se requieren las siguientes herramientas o sus equivalentes:

  * Entorno Unix/Linux.
  * FastQC
  * Fastp
  * Aligner (HISAT2 / STAR)
  * featureCounts
  * R v4.0+ con paquetes: `DESeq2`, `ggplot2`, `pheatmap`.
