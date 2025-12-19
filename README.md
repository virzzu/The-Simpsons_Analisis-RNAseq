# 🧬 Análisis de Expresión Diferencial RNA-seq: Obesidad en personajes de Los Simpsons

Este repositorio contiene el material y los resultados de la actividad práctica sobre **Análisis de Expresión Diferencial (RNA-seq)**. El objetivo del proyecto es identificar genes relacionados con la obesidad comparando perfiles simulados de diferentes personajes de *Los Simpson*.

👉🏼 Puedes ver todo el trabajo redactado y detallado en el siguiente enlace de GitHub Pages: https://virzzu.github.io/The-Simpsons_Analisis-RNAseq/

## 🎯 Objetivos

* Realizar un flujo de trabajo completo de bioinformática (QC, alineamiento, cuantificación).
* Identificar genes diferencialmente expresados (DEGs) utilizando **EdgeR, DESeq2 y limma-voom**.
* Interpretar biológicamente los resultados mediante un **enriquecimiento Funcional**.
* Generar visualizaciones gráficas sobre los resultados obtenidos e interpretarlas.

## 🍩 Datos y Diseño Experimental

Se han propuesto datos simulados de 8 muestras agrupadas en tres fenotipos metabólicos:

| Grupo | Fenotipo | Muestras (Personajes) |
| :--- | :--- | :--- |
| **Obeso 1** | Obesidad Mórbida | Homer, Abraham |
| **Obeso 2** | Sobrepeso/Obesidad | Marge, Patty, Selma |
| **Normopeso** | Control | Bart, Lisa, Maggie |

*Sin embargo, por la organización de la actividad docente, en este repositorio solo se evalúan las muestras correspondientes	a Obeso tipo 2 y Normopeso.

## 📂 Estructura del Repositorio

A continuación se detalla la organización de los archivos:

```text
/
├── 📁 docs/                  # Documentos para generar GitHub pages
│ 
├── 📁 RNAseq/                # Análisis RNAseq
│   ├── Calidad_y_cuantificacion/    # Scripts, datos y resultados para cada paso
│   └── Expresion_diferencial/		del análisis
│   └── Enriquecimiento/
│ 
├── 📁 Poster/                # Poster sobre el análisis y sus documentos
│ 
└── 📄 README.md              # Este archivo
```


Autores: Rita Pellissa Valera, Samuel Pintos González, Tamara Noya Mosquera, Teresa Carrión Mera, Vanesa de las Heras Hermosilla, Virginia García-Loygorri Arias y Yannis Avlonitis Egea.
