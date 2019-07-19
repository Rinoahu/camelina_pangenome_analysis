#!/bin/sh

if [ $# -lt 2 ]

	then
		echo '#######################################'
		echo '#'
		echo '# usage:'
		echo '# this_script.sh [ref.fsa] [qry.list.txt]'
		echo '#'
		echo '#######################################'
		echo ''
		echo '[ref.fsa] : the reference fasta file'
		echo '[qry.list.txt] : the qry.list.txt contain the path of fastq/fasta file the pair-end fastq/fasta should in the same row delim by blank in qry.list.txt [dbsnp]'
		exit 1

fi

# get the name of reference sequence
ref=$1
fasta=`echo $ref | rev | cut -d \/ -f 1 | rev`



ln -sf $fasta $ref.fa
echo "ln -s $fasta $ref.fa"

# create dict for ref
echo "gatk CreateSequenceDictionary -R $ref.fa -O $ref.dict"
gatk CreateSequenceDictionary -R $ref.fa -O $ref.dict

ref=$ref.fa


qry=$2

dbsnp=$3


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
if [ -d "gatk_tmp_dir" ]; then
        echo 'gatk_tmp_dir exist'
        rm -rf gatk_tmp_dir
fi

mkdir gatk_tmp_dir
tmp=gatk_tmp_dir


# creat tmp dir
#for fastq in `cat $qry`
while read fastq

	do

		# 1. bwa mapping and convert to sam step

		# get the reads
		left=`echo $fastq | awk '{print $1}'`
		right=`echo $fastq | awk '{print $2}'`

		# get sample name
	        sample=`echo $left | rev | cut -d \/ -f 1 | rev`
        	sample2=`echo $right | rev | cut -d \/ -f 1 | rev`

		echo 'the left', $left, 'the right', $right

		# mapping the qry to ref and convert sai to sam
		if [ "$right" == "" ]; then
			echo "bwa mem -t 16 -R "@RG\tID:$sample\tLB:$sample\tPL:ILLUMINA\tSM:$sample" $ref $left > $tmp/$sample\_aln.sai.sam"
			bwa mem -t 16 -R "@RG\tID:$sample\tLB:$sample\tPL:ILLUMINA\tSM:$sample" $ref $left > $tmp/$sample\_aln.sai.sam

		else
			echo "bwa mem -t 16 -R "@RG\tID:$sample\tLB:$sample\tPL:ILLUMINA\tSM:$sample" $ref $left $right >  $tmp/$sample\_aln.sai.sam"
			bwa mem -t 16 -R "@RG\tID:$sample\tLB:$sample\tPL:ILLUMINA\tSM:$sample" $ref $left $right >  $tmp/$sample\_aln.sai.sam
        	fi


		# 2. sort aln sam and convert sam to bam
		input=$tmp/$sample\_aln.sai.sam

		echo "gatk SortSam I=$input O=$input.sort SO=coordinate VALIDATION_STRINGENCY=LENIENT"
		gatk SortSam -I $input -O $input.srt -SO coordinate -VALIDATION_STRINGENCY LENIENT


		# 3. delete duplicates
		echo "gatk MarkDuplicates I=$input O=$input.dedup METRICS_FILE=$input.metrics.txt VALIDATION_STRINGENCY=LENIENT"
		gatk MarkDuplicates -I $input.srt -O $input.srt.dedup -METRICS_FILE=$input.srt.metrics.txt -VALIDATION_STRINGENCY=LENIENT



		# 4. add RG to distinguish sample from different source, this step can be avoid if bwa aln has added RG to results
		# index the bam
		echo "samtools index $input.srt.dedup"
		samtools index $input.srt.dedup

		# 5. Indel-based realignment
		if [ "$dbsnp" == "" ]; then

			echo "gatk HaplotypeCaller -R $ref  -I $input -O $input.gatk_raw_gvcf -stand_call_conf 30.0 --emit-ref-confidence GVCF"
			gatk HaplotypeCaller -R $ref -I $input.srt.dedup -O $input.srt.dedup.gatk_raw_gvcf -stand-call-conf 30.0 -ERC GVCF
			echo "gatk GenotypeGVCFs -R $ref --variant $input.srt.dedup.gatk_raw_gvcf -O $input.srt.dedup.gatk_raw_vcf"
			gatk GenotypeGVCFs -R $ref --variant $input.srt.dedup.gatk_raw_gvcf -O $input.srt.dedup.gatk_raw_vcf
			bgzip $input.srt.dedup.gatk_raw_vcf
			tabix -p vcf $input.srt.dedup.gatk_raw_vcf.gz

			echo "samtools mpileup -DSugf $ref $input | bcftools view -Ncvg - > $input.samtools_raw_vcf"
			bcftools mpileup -Sug -f $ref $input.srt.dedup | bcftools call -mv | vcfutils.pl varFilter > $input.srt.dedup.samtools_raw_vcf
			bgzip $input.srt.dedup.samtools_raw_vcf
			tabix -p vcf $input.srt.dedup.samtools_raw_vcf.gz

			echo "vcf-isec $input.srt.dedup.samtools_raw_vcf.gz $input.srt.dedup.gatk_raw_vcf.gz > $input.srt.dedup.vcf"
			vcf-isec $input.srt.dedup.samtools_raw_vcf.gz $input.srt.dedup.gatk_raw_vcf.gz > $input.srt.dedup.vcf


			#echo "gatk SelectVariants -select-type SNP -R $ref --variant $input.srt.dedup.gatk_raw_vcf --concordance $input.srt.dedup.samtools_raw_vcf -O $input.srt.dedup.snp.vcf"
			#gatk SelectVariants -select-type SNP -R $ref --variant $input.srt.dedup.gatk_raw_vcf --concordance $input.srt.dedup.samtools_raw_vcf -O $input.srt.dedup.snp.vcf

			#echo "gatk SelectVariants -select-type INDEL -R $ref --variant $input.srt.dedup.gatk_raw_vcf --concordance $input.srt.dedup.samtools_raw_vcf -O $input.srt.dedup.indel.vcf"
			#gatk SelectVariants -select-type INDEL -R $ref --variant $input.srt.dedup.gatk_raw_vcf --concordance $input.srt.dedup.samtools_raw_vcf -O $input.srt.dedup.indel.vcf

			#echo "gatk MergeVcfs -I $input.srt.dedup.snp.vcf -I $input.srt.dedup.indel.vcf -O $input.srt.dedup.vcf"
			#gatk MergeVcfs -I $input.srt.dedup.snp.vcf -I $input.srt.dedup.indel.vcf -O $input.srt.dedup.vcf

			#echo "vcf"

			# get mean of quality
			MEANQUAL=`awk '{ if ($1 !~ /#/) { total += $6; count++ } } END { print total/count }' $input.srt.dedup.vcf`

			# filter
			echo "gatk VariantFiltration --missing-values-evaluate-as-failing true -V $input.srt.dedup.vcf --verbosity ERROR -R $ref -O $input.srt.dedup.vcf.flt --filter-expression \"QD < 2.0 || FS > 200.0 || SOR > 10.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || QUAL < $MEANQUAL\" --filter-name \"Filter\" "
			gatk VariantFiltration --missing-values-evaluate-as-failing true -V $input.srt.dedup.vcf --verbosity ERROR -R $ref -O $input.srt.dedup.vcf.flt --filter-expression " FS > 200 || SOR > 10 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || QUAL < $MEANQUAL" --filter-name "Filter"


			# add break
			#exit 1

			echo "grep -v Filter $input.srt.dedup.vcf.flt > $input.srt.dedup.vcf.flt.dbsnp"
			grep -v Filter $input.srt.dedup.vcf.flt > $input.srt.dedup.vcf.flt.dbsnp
			dbsnp=$input.srt.dedup.vcf.flt.dbsnp
			gatk IndexFeatureFile -F $dbsnp

		fi


		# 6. Base quality score recalibration
		echo "gatk BaseRecalibrator -R $ref -I $input -O $input.recal_tab --known-sites $dbsnp"
		gatk BaseRecalibrator -R $ref -I $input.srt.dedup -O $input.srt.dedup.recal_tab --known-sites $dbsnp
		gatk ApplyBQSR -R $ref -I $input.srt.dedup --bqsr-recal-file $input.srt.dedup.recal_tab -O $input.srt.dedup.recal_tab.bam

        #echo "gatk PrintReads -R $ref -I $input -BQSR $input.recal_tab -o $input.recal.bam"
		#gatk PrintReads -R $ref -I $input -BQSR $input.recal_tab -o $input.recal_tab.bam

		exit 1

		after_recal_tab=$realn_bam.after_recal_tab
		$gatk_run -T BaseRecalibrator -R $ref -I $input -BQSR $input.recal_tab -o $input.after_recal_tab -knownSites $dbsnp
		echo "$gatk_run -T BaseRecalibrator -R $ref -I $input -BQSR $input.recal_tab -o $input.after_recal_tab -knownSites $dbsnp"

		# 6.1 plot the recalibration
		plots_pdf=$recal_bam\_plots.pdf
		$gatk_run -T AnalyzeCovariates -R $ref -before $input.recal_tab -after $input.after_recal_tab -plots $input.plots.pdf
		echo "$gatk_run -T AnalyzeCovariates -R $ref -before $input.recal_tab -after $input.after_recal_tab -plots $input.plots.pdf"

		# 7. Data compression with reduce reads
		input=$input.recal.bam
		reduced_bam=$recal_bam.reduced_bam
		#$gatk_run -nt 12 -T ReduceReads -R $ref -I $input -o $input.reduced_bam
		#echo "$gatk_run -nt 12 -T ReduceReads -R $ref -I $input -o $input.reduced_bam"

		# 8. Calling variants
		# 8.1 UnifiedGenotyper, UG
		#input=$input.reduced_bam
		ug_output_vcf=$reduced_bam.ug_output_vcf
  		$gatk_run -nt 12 -T UnifiedGenotyper -R $ref -I $input -o $input.ug_output_vcf -stand_call_conf 30 -stand_emit_conf 10
		echo "$gatk_run -nt 12 -T UnifiedGenotyper -R $ref -I $input -o $input.ug_output_vcf -stand_call_conf 30 -stand_emit_conf 10"
		# 8.2 HaplotypeCaller, HC
		hc_output_vcf=$reduced_bam.hc_output_vcf
		$gatk_run -nt 12 -T HaplotypeCall -R $ref -I $input -o $input.hc_output_vcf -stand_call_conf 30 -stand_emit_conf 10 -minPruning 3
		echo "$gatk_run -nt 12 -T HaplotypeCall -R $ref -I $input -o $input.hc_output_vcf -stand_call_conf 30 -stand_emit_conf 10 -minPruning 3"


		# 9. Variant quality score recalibration$input.realn.bam
		input=$input.ug_output_vcf
		recal_file=$ug_output_vcf.recal
		tranches_file=$ug_output_vcf.tranches
		rscript_file=$ug_output_vcf.R
		$gatk_run -nt 12 -T VariantRecalibrator -R $ref -input $input -an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an MQ -mode SNP -recalFile $input.recal -tranchesFile $input.tranches -rscriptFile $input.R -resource:dbsnp,known=true,training=false,truth=false,prior=8.0 $dbsnp

		echo "$gatk_run -nt 12 -T VariantRecalibrator -R $ref -input $input -an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an MQ -mode SNP -recalFile $input.recal -tranchesFile $input.tranches -rscriptFile $input.R -resource:dbsnp,known=true,training=false,truth=false,prior=8.0 $dbsnp"

		recal_vcf=$ug_output_vcf.recal_vcf
		$gatk_run -nt 12 -T ApplyRecalibration -R $ref -input $input -mode SNP -recalFile $input.recal -tranchesFile $input.tranches -o $input.recal_vcf -ts_filter_level 99.0
		echo "$gatk_run -nt 12 -T ApplyRecalibration -R $ref -input $input -mode SNP -recalFile $input.recal -tranchesFile $input.tranches -o $input.recal_vcf -ts_filter_level 99.0"



  		# 10. Genotype Phasing and Reﬁnement
		
		# 11. Functional Annotation

		# 12. Analyzing variant calls

	done < $qry



