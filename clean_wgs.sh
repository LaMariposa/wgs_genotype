#!/bin/bash

#clean_reads.sh
#script to clean raw illumina paired end reads
#M. Supple
#last modified 31 January 2016

#usage clean_reads.sh </path/in/dir>

#requires
	#Trimmomatic
	#FLASh

#input: fastq files (2 for each unit). assumes end in fastq.gz, and R1 and R2 indicate pairs
#output: fastq files--paired (2 for each unit), orphan (1 for each unit), prefix from before R1/R2 in original file names

#get input information
indir=$1

#get file list
file1=($(ls $indir/*R1*fastq.gz))
file2=($(ls $indir/*R2*fastq.gz))


#check if have equal numbers of read1 and read2
if [ ${#file1[@]} = ${#file2[@]} ] 
	then
		echo "processing ${#file1[@]} file pairs"
else 
	echo "Failed!  Number of R1 and R2 files differ"; exit 1
fi


#check for output directory
outdir=clean_fq
if [ ! -d "$outdir" ]; then
  mkdir $outdir
fi


#loop over file pairs and process
for ((a=0; a<${#file1[@]}; a++))
	do
		echo "cleaning file pair ${file1[a]} and ${file2[a]}"

		#determine output prefix
		tempid=`basename ${file1[a]}`
		fileid=${tempid%_R1*}
		echo -e "\n\n$fileid"


		#adapter trim with Trimmomatic
		java -jar /soe/megan/bin/trimmomatic.jar PE -threads 4 ${file1[a]} ${file2[a]} -baseout $outdir/trimmo.fastq.gz \
			ILLUMINACLIP:/soe/megan/bin/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10:1


		#merge overlapping reads with FLASh
		flash -M 140 -t 4 -z -d $outdir $outdir/trimmo_1P.fastq.gz $outdir/trimmo_2P.fastq.gz
		rm $outdir/trimmo_1P.fastq.gz $outdir/trimmo_2P.fastq.gz


		#quality filter and remove small reads with Trimmomatic
		#remaining pairs from flash
		java -jar /soe/megan/bin/trimmomatic.jar PE -threads 4 $outdir/out.notCombined_1.fastq.gz $outdir/out.notCombined_2.fastq.gz \
			-baseout $outdir/trimmo2p.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
		rm $outdir/out.notCombined_1.fastq.gz $outdir/out.notCombined_2.fastq.gz
		#merged from flash
		java -jar /soe/megan/bin/trimmomatic.jar SE -threads 4 $outdir/out.extendedFrags.fastq.gz $outdir/trimmo2u1.fastq.gz \
			LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
		rm $outdir/out.extendedFrags.fastq.gz
		#orphans from trimmomatic
		java -jar /soe/megan/bin/trimmomatic.jar SE -threads 4 $outdir/trimmo_1U.fastq.gz $outdir/trimmo2u2.fastq.gz \
			LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
		java -jar /soe/megan/bin/trimmomatic.jar SE -threads 4 $outdir/trimmo_2U.fastq.gz $outdir/trimmo2u3.fastq.gz \
			LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
		rm $outdir/trimmo_1U.fastq.gz $outdir/trimmo_2U.fastq.gz

		#rename files and merge orphans
		mv $outdir/trimmo2p_1P.fastq.gz $outdir/$fileid\_R1.fastq.gz
                mv $outdir/trimmo2p_2P.fastq.gz $outdir/$fileid\_R2.fastq.gz 		
		cat $outdir/trimmo2p_1U.fastq.gz $outdir/trimmo2p_2U.fastq.gz $outdir/trimmo2u1.fastq.gz $outdir/trimmo2u2.fastq.gz $outdir/trimmo2u3.fastq.gz > $outdir/$fileid\_S.fastq.gz
		rm $outdir/trimmo2p_1U.fastq.gz $outdir/trimmo2p_2U.fastq.gz $outdir/trimmo2u1.fastq.gz $outdir/trimmo2u2.fastq.gz $outdir/trimmo2u3.fastq.gz $outdir/out.histogram $outdir/out.hist

	done

echo -e "\nDONE!!!"


















