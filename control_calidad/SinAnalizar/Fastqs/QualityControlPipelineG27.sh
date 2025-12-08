
#!/bin/bash

SECONDS=0
# Notas: tener instalada una distribucion de conda para poder ejecutar herramientas para el analisis
# Buenas practicas: ejecutar la ayuda de los programas para ver las opciones y sintaxis con programa -h o programa --help
# cambiamos el directorio a donde necesitemos y tengamos nuestros archivos
# cd /Users/virginia/Desktop/Master/SOPG/The-Simpson_Analisis-RNAseq

# # 1. Creamos un entorno en conda con sus programas para utilizarlos de manera segura y aislada de otros entornos
# conda init
# conda create -n Tema4 -c bioconda -c conda-forge -c defaults -c r fastqc fastp multiqc bwa samtools=1.19 \
# unicycler htslib openjdk=17 bandage quast qualimap prokka
# echo "Conda env created\n"

# Sin embargo, ya deberiais tener el entorno, creado. Para activarlo:
# conda env list # aqui veis las listas de entornos que teneis dentro de conda
# conda activate Tema4 # lo mas probable, segun lo que hemos dado en clase, sea que se llame asi

# Con este comando el entorno Tema4 y le instalamos los contenedores, que son gestores de paquetes (-c) bioconda, conda-forge, defaults y r. 
# Posteriormente le mete los programas con sus versiones:
# - fastqc para el control de calidad
# - fastp para filtrar por longitud y calidad de secuencias
# - multiqc para informes de fastqc multiples
# - bwa para mappear las secuencias
# - samtools para trabajar con ficheros del mappeo
# - unicycler para ensamblado de organismos pequenios
# - htslib y openjdk son dependencias de las anteriores
# - bandage para visualizar los ensamblados
# - quast para control de calidad de ensamblado
# - qualimap para calidad de mapeos
# - prokka para anotacion/notacion de genes

# 2. Activamos el entorno que acabamos de crear, y automaticamente nos metemos dentro de el
# conda activate Tema4a
# echo "Conda env activated!\n"

# nos movemos a la carpeta contenedora de los fastq
cd ./Fastqs

# creamos las carpetas que vamos a necesitar
mkdir -p Quality/Raw Quality/Filtered Trimmed

# En qualityraw haremos los controles de calidad de cada uno de los crudos con fastqc
# y luego cada uno lo meteremos en los trimmed con lo que devuelva fastp,
# finalmente meteremos en qualityfiltered lo que devuelva fastqc sobre los archivos de trimmed

###########################################################################################
###########################################################################################
############                                                                   ############
######                                                                               ######
###### !!! OJO, borramos todos los archivos que no nos corresponden (homer, abraham) ######
######                                                                               ######
############                                                                   ############
###########################################################################################
###########################################################################################


# 3. Realizamos el control de calidad de todos los archivos con el programa FASTQC ( -t es el n de hilos que ejecutara el prog)
fastqc *fastq.gz -o Quality/Raw/ -t 32

echo "Quality control done by fastqc, outputted in Quality/Raw\n"

# 4. Ejecutamos el filtrado de secuencias que queremos por calidad y longitud con FASTP para cada muestra con un bucle for
ls *fastq.gz | cut -d _ -f 1 | sort -u > muestras.txt
# aqui metemos los nombres de las muestras en una lista para iterar sobre ello despues

for i in $(cat muestras.txt); do fastp --in1 $i*R1* --in2 $i*R2* \
--out1 Trimmed/$i"_R1_filtered.fastq.gz" --out2 Trimmed/$i"_R2_filtered.fastq.gz" \
--detect_adapter_for_pe --cut_front --cut_tail --cut_window_size 12 --cut_mean_quality 30 --length_required 35 \
--json Trimmed/$i.json --html Trimmed/$i.html --thread 32; done

echo "Filtering done by fastp. Output can be found in Trimmed/ \n"

# 5. Ejecutamos de nuevo un segundo control de calidad para corroborar que hemos realizado correctamente el filtrado
fastqc Trimmed/*fastq.gz -o Quality/Filtered/ --threads 32

echo "Second quality control to filter sequences from fastp done by fastqc. Outputted in Quality/Filtered/\n"

duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed.\n"

# tambien podemos generar un informe resumen con los pasos que hemos hecho con multiqc, se aplica directamente a la carpeta contenedora de todos los informes y el html lo deja en la misma
multiqc .

open multiqc_report.html

