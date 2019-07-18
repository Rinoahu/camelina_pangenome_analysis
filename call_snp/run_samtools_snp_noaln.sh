#!/bin/sh

if [ $# -lt 2 ]

	then
		echo 'this script is used to call snp from bam files'
		echo 'usage:'
		echo 'this_script.sh ref.fsa bam.list.txt'
		echo 'the qry.list.txt contain the path of bam file'
		exit 1

fi

ref=$1
qry=$2


if [ -d 'samtools_tmp_dir' ];
	then
		rm -rf samtools_tmp_dir
fi

mkdir samtools_tmp_dir

# index the ref
samtools faidx $ref


# creat tmp dir
for bam in `cat $qry`

	do

		echo 'bam file', $bam
		left=$bam

		# sort the snp
		samtools sort -@ 16 $left $left.sort

		# reduce dup
		samtools rmdup $left.sort.bam $left.sort.bam.rmdup_bam

		total_bams=`echo  $total_bams $left.sort.bam.rmdup_bam`
	done


samtools mpileup -m 3 -F 0.0002 -uf $ref $total_bams | bcftools view -bvcg - > ./samtools_tmp_dir/total_bams.raw.bcf
bcftools view ./samtools_tmp_dir/total_bams.raw.bcf > ./samtools_tmp_dir/total_bams.raw.bcf.vcf
bcftools view ./samtools_tmp_dir/total_bams.raw.bcf | vcfutils.pl varFilter -w 3 -d 10 -Q 20 > ./samtools_tmp_dir/total_bams.raw.bcf.fl.vcf

# clean the temp files
mv *_aln.sai* samtools_tmp_dir

