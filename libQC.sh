#!/bin/bash

#libQC.sh
#script to QC MiSeq sequenced libraries
#M. Supple
#last modified 9 Nov 2017

#usage align_reads.sh </path/genome/genome.fasta> </path/in/dir>

#requires
	#fastqc
	#bwa

#input: paired fastq files (assumes *_R1_001.fastq.gz and *_R2_001.fastq.gz)
#output:
	#fastqc report 
	#sam file and fragment length distribution
	#prinseq output
	#flagstat

threads=6

#get input information
genome=$1
indir=$2

#get file list
file1=($(ls $indir/*_R1_001.fastq.gz))
file2=($(ls $indir/*_R2_001.fastq.gz))

#array of library names
libsids=()

#check if have equal numbers of read1 and read2
if [ ${#file1[@]} = ${#file2[@]} ] 
	then
		echo -e "processing ${#file1[@]} file pairs\naligning to genome $genome"
else 
	echo "Failed!  Number of R1 and R2 files differ"; exit 1
fi


#check for output directory
fqcdir=fastqc
if [ ! -d "$fqcdir" ]; then
  mkdir $fqcdir
fi

aligndir=align_bam
if [ ! -d "$aligndir" ]; then
  mkdir $aligndir
fi

tlendir=tlens
if [ ! -d "$tlendir" ]; then
  mkdir $tlendir
fi

psdir=prinseq
if [ ! -d "$psdir" ]; then
  mkdir $psdir
fi


#loop over file pairs and process
for ((a=0; a<${#file1[@]}; a++))
	do
		echo -e "\n\nQCing file pair ${file1[a]} and ${file2[a]}"

		#fastqc
		fastqc -o $fqcdir ${file1[a]}
		fastqc -o $fqcdir ${file2[a]}	

		#determine output prefix
		tempid=`basename ${file1[a]}`
		fileid=${tempid%_R1*}
		echo "read group=$fileid"

		#align with bwa mem
		bwa mem -t $threads -PM $genome ${file1[a]} ${file2[a]} > $aligndir/$fileid.sam
		echo "\n\n\nstats for $fileid" >> flagstat.out
		samtools flagstat $aligndir/$fileid.sam >> flagstat.out
		
		#get frag length distribution
		cut -f9 -d$'\t' $aligndir/$fileid.sam > $tlendir/$fileid.tlen.txt

		#prinseq
		perl /projects/redser3-notbackedup/projects/pheintzman/Scripts/prinseq-lite.pl -fastq ${file1[a]} -graph_data $psdir/$fileid.R1.gd \
			-graph_stats ld,gc,qd,ns,pt,ts,aq,de,da,sc -out_good null -out_bad null
		perl /projects/redser3-notbackedup/projects/pheintzman/Scripts/prinseq-lite.pl -fastq ${file2[a]} -graph_data $psdir/$fileid.R2.gd \
			-graph_stats ld,gc,qd,ns,pt,ts,aq,de,da,sc -out_good null -out_bad null
		perl /projects/redser3-notbackedup/projects/pheintzman/Scripts/prinseq-graphs.pl -i $psdir/$fileid.R1.gd -o $psdir/$fileid.R1.graph_data -png_all -html_all 
		perl /projects/redser3-notbackedup/projects/pheintzman/Scripts/prinseq-graphs.pl -i $psdir/$fileid.R2.gd -o $psdir/$fileid.R2.graph_data -png_all -html_all
	done

Rscript /soe/megan/wgs_genotype/tlendist.r

echo -e "\n\nDONE!!!"





















