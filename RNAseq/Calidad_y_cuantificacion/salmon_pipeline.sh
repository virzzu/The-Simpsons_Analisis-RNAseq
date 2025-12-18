#!/bin/bash

# ------------------------ Pipeline de pseudoalineamiento y cuantificacion con Salmon ----------------------------
# ----------------------------------------------------------------------------------------------------------------
 # 										!!!!!!  ATENCION  !!!!!!
# Lo primero que hay que hacer es * crear un nuevo entorno * porque con el entorno del Tema 4 hay incompatibilidades y el resultado
# sale mal.

# conda create --name salmon_env 
# conda activate salmon_env
# conda install -c bioconda salmon

# Utilizamos Salmon porque lleva a cabo un pseudoalineamiento, realmente no nos importa tanto la posicion exacta base por base, sino 
# solo saber si la lectura es compatible con alguno de los transcritos de referencia que le pasemos
# Ademas, Salmon no genera archivos .sam o .bam tipicos del mapeo, que son muy pesados, sino que devuelve directamente
# la matriz de conteos de las lecturas.

############# Workflow Salmon #############
# despues de instalar salmon, el workflow de alineamiento y cuantificacion seria:
		# 1- indexar la referencia fasta
		# 2- quantificar los fastq contra el index

# 1- Indexar la referencia fasta
salmon index -t Referencia.fasta -i salmon_index
		# aqui utilizamos la herramienta index, pasandole el fasta que queremos que indexe con -t (de transcrito)
		# y tenemos que indicarle con -i el repo de salida en el que va a sacarlo

# Salmon al ejecutar esto, crea una carpeta salmon_index con varios archivos de funcionamiento interno entre los que se encuentra
# informacion sobre la secuencia, una tabla hash de k-meros, un array de sufijos para buscar subcadenas durante el mapeo,
# y metadatos sobre la secuencia y sobre el analisis (pej version del algoritmo). 


# 2- Cuantificar los fastq frente al index que ha creado salmon
mkdir quant_results
# creamos una carpeta donde guardar la cuantificacion

# con un bucle, quantificamos todas las muestras
for i in $(cat Fastqs/muestras.txt);
do 
	echo "Quantifying sample: ${i}"
	salmon quant -i salmon_index -l A \
		-1 Fastqs/${i}_R1.fastq.gz \
		-2 Fastqs/${i}_R2.fastq.gz \
		--validateMappings \
		-o quant_results/${i}_quant
	echo "Quantifying done"
done
	# -1 pide la primera lectura del paired-end
	# -2 pide la segunda lectura del paired-end
	# --validateMappings es una version mejorada del algoritmo interno
	# -o donde y como queremos que salga el output

# Salmon quant devolvera mensajes que tenemos que leer e interpretar, podemos buscar en la terminal con ctrl+f:
	# [info] mapping rate = % -> aqui deberiamos estar viendo porcentajes del orden de >70-80%. 
	# 	(en el trabajo no sale, puede que sea porque son datos artificiales)
	# [info] counted n total reads -> sirve para saber si tenemos datos suficientes para hacer la estadistica (la quantf en si)
	# [error] o [warning] -> siempre es importante leer

# Asi mismo, los archivos importantes que nos vamos a llevar a R para analizar su diff expresion van a ser los quant.sf

# finalmente podemos lanzar un multiqc con los archivos .sf para que lo explique todo un poco en html
# pero para ello tenemos que cambiar de entorno porque ya no tenemos instalado multiqc aqui :(
