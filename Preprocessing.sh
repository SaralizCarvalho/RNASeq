#!/bin/bash

export PATH=$PATH:~/bin/FastQC
trim_path="$HOME/bin/Trimmomatic-0.39/trimmomatic-0.39.jar"

mkdir -p preprocessing/hq_reads
#mkdir -p preprocessing/fastqc/raw
mkdir -p preprocessing/fastqc/filtered

R1_suf=$2
R2_suf=$3

#first run of fastqc for raw reads

while read f
do
	in1=${f}${R1_suf}
	in2=${f}${R2_suf}
	date
	sample=$(basename $f)
	echo "Running fastqc for sample ${sample}"
	echo "fastqc -t 10 -o preprocessing/fastqc/raw ${in1} ${in2}"
	#fastqc -t 10 -o preprocessing/fastqc/raw ${in1} ${in2}
	
	#Running Trimmomatic
	echo "Running Trimmomatic for sample ${sample}"
	echo "java -jar $trim_path PE $in1 $in2 -baseout preprocessing/hq_reads/${sample}.filtered.fq.gz -threads 20 -trimlog preprocessing/hq_reads/${sample}.filtered.log LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50"
	java -jar $trim_path PE $in1 $in2 -baseout preprocessing/hq_reads/${sample}.filtered.fq.gz -threads 20 -trimlog preprocessing/hq_reads/${sample}.filtered.log LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50

done < $1

list=""
for f in $(ls -1 preprocessing/hq_reads/*P.fq.gz) 
do
	list+=" ${f}"
done

echo "Running FASTQC for filtered reads"
fastqc -t 10 -o preprocessing/fastqc/filtered $list
