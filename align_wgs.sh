#!/bin/bash

#align_reads.sh
#script to align cleand illumina reads
#M. Supple
#last modified 3 February 2016

#usage align_reads.sh </path/in/dir>

#requires
	#bwa
	#samtools
	#picard

#input: fastq files (3 for each unit). assumes end in clean.fastq.gz, and _R1_ and _R2_ indicate pairs, _S_ indicates orphans
#output: bam file and index for each unit

#get input information
indir=$1

#get file list
file1=($(ls $indir/*_R1_clean.fastq.gz))
file2=($(ls $indir/*_R2_clean.fastq.gz))
file3=($(ls $indir/*_S_clean.fastq.gz))


#check if have equal numbers of read1 and read2 and read3
if [ ${#file1[@]} = ${#file2[@]} ] && [ ${#file2[@]} = ${#file3[@]} ] 
	then
		echo "processing ${#file1[@]} file triples"
else 
	echo "Failed!  Number of R1, R2, S files differ"; exit 1
fi


#check for output directory
outdir=align_bam
if [ ! -d "$outdir" ]; then
  mkdir $outdir
fi


#loop over file pairs and process
for ((a=0; a<${#file1[@]}; a++))
	do
		echo -e "\n\ncleaning file triple ${file1[a]} and ${file2[a]} and ${file3[a]}"

		#determine output prefix
		tempid=`basename ${file1[a]}`
		fileid=${tempid%_R1*}
		echo "$fileid"


		#align with bwa mem
		#pairs
		bwa mem -t 6 -PM -R '@RG\tID:test\tSM:test'
		#orphans

		#mark duplicates with picard
		
		#base recalibration
		#realign around indels?

		#index



	done

echo -e "\nDONE!!!"


















