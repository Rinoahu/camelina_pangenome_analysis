#!/bin/sh
###########################################################################################
#
#    Copyright (C) 2014, Hu Xiao <huxiao@sibs.ac.cn>
#
#    The run_samtools_snp is an automate pipeline to call snp by using bwa and samtools
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
############################################################################################


# this shell script is an automate pipeline to call snp by bwa and samtools
#
# usage:
# sh this_script.sh ref.fsa qry.list.txt
# 
#      ref.fsa: the fasta file of genome
# qry.list.txt: the qry.list.txt contain the path of fastq/fasta file
#               the pair-end fastq/fasta should in the same row delim by blank in qry.list.txt'
#


if [ $# -lt 2 ]

	then
		echo 'usage:'
		echo 'sh this_script.sh ref.fsa qry.list.txt'
		echo ''
		echo '     ref.fsa: the fasta file of genome'
		echo 'qry.list.txt: the qry.list.txt contain the path of fastq/fasta file'
		echo '              the pair-end fastq/fasta should in the same row delim by blank in qry.list.txt'
		echo ''
		exit 1

fi

ref=$1
qry=$2


# index the ref
if [ -f "$ref.bwt" ]; then
	echo 'ref has been indexed by bwa'

else
	bwa index $ref

fi

if [ -f "$qry.fai" ]; then

	echo 'fasta of ref has been index by samtools faidx'

else
	samtools faidx $ref

fi


# create tmp dir
if [ -d "samtools_tmp_dir" ]; then
	echo 'samtools_tmp_dir exist'
	rm -rf samtools_tmp_dir
fi

mkdir samtools_tmp_dir

tmp=samtools_tmp_dir

total_bams=""

while read fastq

	do
		left=`echo $fastq | awk '{print $1}'`
		right=`echo $fastq | awk '{print $2}'`

		# get sample name
		sample=`echo $left | rev | cut -d \/ -f 1 | rev`
		sample2=`echo $right | rev | cut -d \/ -f 1 | rev`

		echo 'the left', $left, 'the right', $right

		# mapping the qry to ref and convert sai to sam
		bwa aln -t 16 $ref $left > $tmp/$sample\_aln.sai
		if [ "$right" == "" ]; then
			bwa samse -r "@RG\tID:$sample\tLB:$sample\tPL:ILLUMINA\tSM:$sample" $ref $tmp/$sample\_aln.sai $left  > $tmp/$sample\_aln.sai.sam

		else
			bwa aln -t 16 $ref $right > $tmp/$sample2\_aln.sai
			bwa sampe -r "@RG\tID:$sample\tLB:$sample\tPL:ILLUMINA\tSM:$sample" $ref $tmp/$sample\_aln.sai $tmp/$sample2\_aln.sai $left $right >  $tmp/$sample\_aln.sai.sam

		fi

		# convert sam to bam
		samtools view -bS -@ 16 $tmp/$sample\_aln.sai.sam > $tmp/$sample\_aln.sai.sam.bam

		# sort the snp
		samtools sort -m 40G -@ 16 $tmp/$sample\_aln.sai.sam.bam $tmp/$sample\_aln.sai.sam.bam.sort

		# reduce dup
		if [ "$right" == "" ]; then
			samtools rmdup -s $tmp/$sample\_aln.sai.sam.bam.sort.bam $tmp/$sample\_aln.sai.sam.bam.sort.bam.rmdup_bam
		else
			samtools rmdup -S $tmp/$sample\_aln.sai.sam.bam.sort.bam $tmp/$sample\_aln.sai.sam.bam.sort.bam.rmdup_bam
		fi

		total_bams=`echo $total_bams $tmp/$sample\_aln.sai.sam.bam.sort.bam.rmdup_bam`

	done < $qry

echo 'call snp', $bams
echo 'the ref', $ref

samtools mpileup -uf $ref $total_bams | bcftools view -bvcg - > $tmp/total_bams.raw.bcf
bcftools view $tmp/total_bams.raw.bcf | vcfutils.pl varFilter -w 5 -d 10 -Q 20 > $tmp/total_bams.raw.bcf.fl.vcf
awk '$6 >= 20' $tmp/total_bams.raw.bcf.fl.vcf > $tmp/total_bams.raw.bcf.fl.vcf.qual20

if [ -d "samtools_results_dir" ]; then
	rm -rf samtools_results_dir

fi

mv samtools_tmp_dir samtools_results_dir
# clean the temp files
#mkdir tmp_dir
#mv *_aln.sai* tmp_dir

