#!/bin/bash

# This script does the mapping of paired end reads to a reference genome using STAR

export PATH=$PATH:$HOME/bin/STAR-2.7.9a/bin/Linux_x86_64
export PATH=$PATH:$HOME/bin/subread-2.0.3-Linux-x86_64/bin

# Specify directories
gindex="$HOME/data/sara/genome_index"
gannotation="$HOME/data/sara/genome/GCF_002906115.3_Cork_oak_2.0_genomic.gtf"
genome="$HOME/data/sara/genome/GCF_002906115.3_Cork_oak_2.0_genomic.fna"
readdir="$HOME/data/sara/corkoak_reads/preprocessing/hq_reads/"
outdir="$HOME/data/sara/corkoak_reads/aligned_reads"

mkdir -p $outdir
mkdir -p $gindex


# Specify suffix
r1_suf=".filtered_1P.fq.gz"
r2_suf=".filtered_2P.fq.gz"


# Running STAR to create genome index
# Only run this line the FIRST time for your data

#STAR --runMode genomeGenerate --genomeDir $gindex --genomeFastaFiles $genome --runThreadN 20 --sjdbGTFfile $gannotation --genomeSAindexNbases 13

# Running STAR for sample list (one sample per line)

while read -r sample
do
	in1="${readdir}${sample}${r1_suf}"
        echo "$in1"
	in2="${readdir}${sample}${r2_suf}"
	echo "$in2"

	echo "Running STAR for sample ${sample}"
	echo "STAR --runThreadN 20 --genomeDir $gindex --readFilesIn $in1 $in2 --outFileNamePrefix ${outdir}/${sample} --twopassMode Basic --readFilesCommand gunzip -c --outReadsUnmapped Fastx --outFilterIntronMotifs RemoveNoncanonical --alignIntronMax 100000 --outSAMtype BAM Unsorted"

	STAR --runThreadN 20 --genomeDir $gindex --readFilesIn $in1 $in2 --outFileNamePrefix ${outdir}/${sample} --twopassMode Basic --readFilesCommand gunzip -c --outReadsUnmapped Fastx --outFilterIntronMotifs RemoveNoncanonical --alignIntronMax 100000 --outSAMtype BAM SortedByCoordinate  

done < $1

#Prepare BAM list

outbam="$outdir"/"*Aligned.out.bam"
bamlist=""
for bam in $(ls -1 $outbam); do bamlist+=" ${bam}"; done

# Running FeatureCounts to obtain a gene count table

echo "featureCounts -T 10 -p --countReadPairs -t exon -g gene -a $gannotation -o ${outdir}/raw_counts.txt $bamlist"

featureCounts -T 10 -p --countReadPairs -t exon -g gene -a $gannotation -o ${outdir}/raw_count$bamlist


