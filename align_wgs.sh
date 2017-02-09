#!/bin/bash

#align_reads.sh
#script to align cleand illumina reads
#M. Supple
#last modified 3 February 2016

#usage align_reads.sh </path/genome/genome.fasta> </path/in/dir>

#requires
	#bwa
	#samtools
	#picard

#input: fastq files (3 for each unit). assumes end in clean.fastq.gz, and _R1_ and _R2_ indicate pairs, _S_ indicates orphans
#output: bam file and index for each library

threads=6
picard=/soe/megan/bin/picard.jar
gatk=/soe/megan/bin/GenomeAnalysisTK.jar 

#get input information
genome=$1
indir=$2

#get file list
file1=($(ls $indir/*_R1_clean.fastq.gz))
file2=($(ls $indir/*_R2_clean.fastq.gz))
file3=($(ls $indir/*_S_clean.fastq.gz))

#array of library names
libsids=()

#check if have equal numbers of read1 and read2 and read3
if [ ${#file1[@]} = ${#file2[@]} ] && [ ${#file2[@]} = ${#file3[@]} ] 
	then
		echo -e "processing ${#file1[@]} file triples\naligning to genome $genome"
else 
	echo "Failed!  Number of R1, R2, S files differ"; exit 1
fi


#check for output directory
outdir=align_bam
if [ ! -d "$outdir" ]; then
  mkdir $outdir
fi


#loop over file triples and process
for ((a=0; a<${#file1[@]}; a++))
	do
		echo -e "\n\ncleaning file triple ${file1[a]} and ${file2[a]} and ${file3[a]}"

		#determine output prefix
		tempid=`basename ${file1[a]}`
		fileid=${tempid%_R1*}
		echo "read group=$fileid"
		#determine library id
		libid=${fileid%_S*}
		echo "library ID=$libid"
		libsids+=($libid)

		#align with bwa mem
		#pairs
		bwa mem -t $threads -PM -R "@RG\tID:$fileid\tLB:$libid\tSM:$libid\tPL:ILLUMINA" $genome \
			${file1[a]} ${file2[a]} > $outdir/temp_p.sam
		#orphans
		bwa mem -t $threads -M -R "@RG\tID:$fileid\tLB:$libid\tSM:$libid\tPL:ILLUMINA" $genome \
                        ${file3[a]} > $outdir/temp_u.sam

		#generate bams and merge pairs and orphans and sort
		samtools view -@ $threads -b -o $outdir/temp_p.bam $outdir/temp_p.sam
		samtools view -@ $threads -b -o $outdir/temp_u.bam $outdir/temp_u.sam
		samtools merge -c $outdir/temp.bam $outdir/temp_p.bam $outdir/temp_u.bam
		samtools sort -@ $threads -o $outdir/$fileid.bam -T tmpsort $outdir/temp.bam

		#clean
		rm $outdir/temp_p.sam $outdir/temp_u.sam $outdir/temp_p.bam $outdir/temp_u.bam $outdir/temp.bam

	done



#loop over libraries
uniq_libsids=($(printf "%s\n" "${libsids[@]}" | sort -u))
echo -e "\n\n${#uniq_libsids[@]} unique libraries: ${uniq_libsids[@]}"
for ((b=0; b<${#uniq_libsids[@]}; b++))
	do
		echo -e "\n\nprocessing library ${uniq_libsids[b]}"

		#merge bams from each library and mark duplicates
		#get a list of all the bams for the library
		libbams=`echo $outdir/${uniq_libsids[b]}*.bam | sed 's/ / I=/g'`
		echo "using bam files $libbams"
		java -jar $picard MarkDuplicates M=$outdir/${uniq_libsids[b]}.metrics I=$libbams O=$outdir/${uniq_libsids[b]}_markdup.bam
	
		#index
		samtools index $outdir/${uniq_libsids[b]}_markdup.bam

		#call library variants
		java -jar $gatk -T HaplotypeCaller -R $genome -I $outdir/$fileid.bam --emitRefConfidence GVCF -o $outdir/$fileid.g.vcf

	done


echo -e "\n\nDONE!!!"


















