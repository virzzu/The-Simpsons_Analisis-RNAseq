---
title: "Inicio"
output: 
  html_document: 
    toc: false
    highlight: tango
    css: style.css
---
Este proyecto presenta un flujo de trabajo bioinformático completo para el análisis de expresión diferencial de genes (RNA-seq). El objetivo del estudio es identificar firmas transcriptómicas y marcadores moleculares alterados al comparar sujetos con Normopeso frente a sujetos con Sobrepeso Tipo 2.

Se ha implementado un pipeline en lenguajes bash y R que abarca desde el preprocesamiento de las librerías y generación de la matriz de conteos, la normalización de los conteos, la visualización avanzada de resultados y su identificación y significado biológico. Mediante el uso de modelos estadísticos (EdgeR), se detectaron genes clave asociados a la regulación metabólica (como FTO y LEP, entre otros), permitiendo caracterizar el perfil molecular distintivo del sobrepeso.

❗️ Este análisis se ha realizado a partir de un dataset de muestras simuladas generadas con fines puramente académicos y didácticos. Sin embargo, los resultados muestran coherencia biológica con la literatura científica sobre obesidad y metabolismo. El propósito principal de este repositorio es demostrar la competencia técnica en el manejo de flujos de trabajo de análisis transcriptómico.

### Objetivos de aprendizaje
 * Realizar un análisis de calidad y alineamiento de lecturas RNA-seq.
 * Cuantificar niveles de expresión génica y normalizar los datos.
 * Identificar genes diferencialmente expresados entre grupos de individuos.
 * Interpretar perfiles de expresión en función del fenotipo o condición.
 * Representar los resultados en un póster científico.
