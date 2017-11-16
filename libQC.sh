#!/bin/bash

#libQC.sh
#script to QC MiSeq sequenced libraries
#M. Supple
#last modified 9 Nov 2017

#usage align_reads.sh </path/genome/genome.fasta> </path/in/dir>

#requires
	#FastQC
	#BWA
	#SAMtools 
	#PRINSEQ
	#Trimmomatic
	#FLASh

#input: paired fastq files (assumes *_R1_001.fastq.gz and *_R2_001.fastq.gz)

#output:
	#fastqc reports 
	#fragment length distribution from alignment
	#flagstat output
	#prinseq reorts
	#stats from adapter trimming and merging

threads=6
adapters=/soe/megan/bin/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa
trimmo=/soe/megan/bin/trimmomatic.jar

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

tlendir=tlens
if [ ! -d "$tlendir" ]; then
  mkdir $tlendir
fi

psdir=prinseq
if [ ! -d "$psdir" ]; then
  mkdir $psdir
fi

tempdir=temp
if [ ! -d "$tempdir" ]; then
  mkdir $tempdir
fi

#loop over file pairs and process
for ((a=0; a<${#file1[@]}; a++))
	do
		echo -e "\n\nQCing file pair ${file1[a]} and ${file2[a]}"

		#fastqc
		fastqc -dir $tempdir -o $fqcdir ${file1[a]}
		fastqc -dir $tempdir -o $fqcdir ${file2[a]}	

		#determine output prefix
		tempid=`basename ${file1[a]}`
		fileid=${tempid%_R1*}
		echo "read group=$fileid"

		#align with bwa mem
		bwa mem -t $threads -PM $genome ${file1[a]} ${file2[a]} > temp.sam
		echo -e "\n\n\nstats for $fileid" >> flagstat.out
		samtools flagstat temp.sam >> flagstat.out
		
		#get fragment length distribution
		cut -f9 -d$'\t' temp.sam > $tlendir/$fileid.tlen.txt
		rm temp.sam

		#prinseq
		zcat ${file1[a]} | prinseq-lite.pl -fastq stdin -graph_data temp.R1.gd \
			-graph_stats ld,gc,qd,ns,pt,ts,aq,de,da,sc -out_good null -out_bad null
		zcat ${file2[a]} | prinseq-lite.pl -fastq stdin -graph_data temp.R2.gd \
                       -graph_stats ld,gc,qd,ns,pt,ts,aq,de,da,sc -out_good null -out_bad null
		prinseq-graphs.pl -i temp.R1.gd -o $psdir/$fileid.R1.graph_data -html_all 
		prinseq-graphs.pl -i temp.R2.gd -o $psdir/$fileid.R2.graph_data -html_all

		#adapter trim with trimmomatic
		java -jar $trimmo PE -threads 4 ${file1[a]} ${file2[a]} -baseout temp \
			ILLUMINACLIP:$adapters:2:30:10:1

		#merge reads with flash
		flash -M 140 -t 4 -z temp_1P temp_2P

		#cleanup
		rm temp.R1.gd temp.R2.gd
		rm temp_1P temp_2P temp_1U temp_2U
		rm out.notCombined_1.fastq.gz out.notCombined_2.fastq.gz out.extendedFrags.fastq.gz 
		rm out.hist out.histogram
		
	done

#make fragment length distributino
Rscript /soe/megan/wgs_genotype/tlendist.r
#get Trimmomatic stats
grep 'read group=\|Input Read Pairs:' qc.out > cleaning.stats
#get FLASh stats
grep 'read group=\|Percent combined:' qc.out >> cleaning.stats

#more cleaning
rm -I fastqc/*zip
rmdir temp

echo -e "\n\nDONE!!!"





















