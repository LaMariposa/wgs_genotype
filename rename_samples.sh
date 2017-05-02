#!/bin/bash

#rename_samples.sh
#script to change read group sample names
#M. Supple
#last modified 17 February 2017

#usage ./rename_samples.sh <input.txt> <bam/dir>

#requires
	#samtools

#input: bam files and file containing current names (seqid) and sample names (sampid)
#output: bam file and index for each sample

samtools=/soe/megan/bin/samtools

#read input info
source $1
indir=$2

#check if have equal numbers of seqid and sampid
if [ ${#seqid[@]} = ${#sampid[@]} ]
	then
		echo -e "processing ${#seqid[@]} file pairs"
else 
	echo "Failed!  Different number of current and new names"; exit 1
fi


#check for output directory
outdir=bam_rename
if [ ! -d "$outdir" ]; then
  mkdir $outdir
fi


#loop over files and process
for ((a=0; a<${#seqid[@]}; a++))
	do

		echo -e "\n\nrenaming ${seqid[a]} to ${sampid[a]}"
		
		$samtools view -H $indir/${seqid[a]}_markdup.bam | sed "s/SM:${seqid[a]}/SM:${sampid[a]}/g" \
			| $samtools reheader - $indir/${seqid[a]}_markdup.bam > $outdir/${sampid[a]}_${seqid[a]}_markdup.bam

	done





echo -e "\n\nDONE!!!"


















