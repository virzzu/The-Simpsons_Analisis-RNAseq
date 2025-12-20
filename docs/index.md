---
title: "Introducción"
output: 
  html_document: 
    toc: false
    highlight: tango
    css: style.css
---
##

La obesidad es un trastorno metabólico complejo influenciado por múltiples factores genéticos y ambientales. Para comprender mejor los mecanismos moleculares de esta condición, en este proyecto se ha llevado a cabo un análisis de expresión diferencial utilizando datos de secuenciación de RNA (RNA-seq).  Para el estudio, se ha trabajado con un conjunto de datos simulados pertenecientes a 8 personajes del universo de Los Simpson, los cuales han sido agrupados en tres perfiles según su condición metabólica: Grupo Obeso 1 (Abraham y Homer Simpson), Grupo Obeso 2 (Marge, Patty y Selma Bouvier) y Grupo Normopeso (Bart, Lisa y Maggie Simpson). Esta configuración permite modelar diferencias genéticas y evaluar cambios en la expresión.

La metodología aplicada abarca el control de calidad de los archivos FASTQ, el alineamiento y cuantificación de lecturas, así como la necesaria normalización de los datos para eliminar sesgos técnicos. Posteriormente, se ha realizado el análisis estadístico de expresión diferencial para detectar genes sobreexpresados o infraexpresados en las comparativas asignadas. Finalmente, para aportar un contexto funcional, se ha incluido un análisis de enriquecimiento de rutas metabólicas y procesos biológicos, buscando relacionar los genes alterados con mecanismos fisiopatológicos conocidos de la obesidad.

Este proyecto presenta un flujo de trabajo bioinformático completo para el análisis de expresión diferencial de genes (RNA-seq). El objetivo del estudio es identificar firmas transcriptómicas y marcadores moleculares alterados al comparar sujetos con Normopeso frente a sujetos con Sobrepeso Tipo 2.

Se ha implementado un pipeline en lenguajes bash y R que abarca desde el preprocesamiento de las librerías y generación de la matriz de conteos, la normalización de los conteos, la visualización avanzada de resultados y su identificación y significado biológico. Mediante el uso de modelos estadísticos (EdgeR), se detectaron genes clave asociados a la regulación metabólica (como FTO y LEP, entre otros), permitiendo caracterizar el perfil molecular distintivo del sobrepeso tipo 2.

❗️ Este análisis se ha realizado a partir de un dataset de muestras simuladas generadas con fines puramente académicos y didácticos. Sin embargo, los resultados muestran coherencia biológica con la literatura científica sobre obesidad y metabolismo. El propósito principal de este repositorio es demostrar la competencia técnica en el manejo de flujos de trabajo de análisis transcriptómico.

**Objetivos de aprendizaje**

 * Realizar un análisis de calidad y alineamiento de lecturas RNA-seq.
 
 * Cuantificar niveles de expresión génica y normalizar los datos.
 
 * Identificar genes diferencialmente expresados entre grupos de individuos.
 
 * Interpretar perfiles de expresión en función del fenotipo o condición.
 
 * Representar los resultados en un póster científico.
