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

# set the java home and the picard home
JAVA_HOME=''

picard_run='java -Xmx4g -jar /home/zhans/tools/assembly/broad/picard/dist'
gatk_run='java -Xmx4g -jar /home/zhans/tools/assembly/broad/gatk/GenomeAnalysisTK.jar'




ref=$1
fasta=`echo $ref | rev | cut -d \/ -f 1 | rev`
# create dict for ref
$picard_run/CreateSequenceDictionary.jar R=$ref O=$ref.dict
echo "$picard_run/CreateSequenceDictionary.jar R=$ref O=$ref.dict"
ln -s $fasta $ref.fa
echo "ln -s $fasta $ref.fa"
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


#$gatk_run/CreateSequenceDictionary.jar I=$ref O=$ref_nosuf.dict


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
		bwa aln -t 16 $ref $left > $tmp/$sample\_aln.sai
		if [ "$right" == "" ]; then
			bwa samse -r "@RG\tID:$sample\tLB:$sample\tPL:ILLUMINA\tSM:$sample" $ref $tmp/$sample\_aln.sai $left > $tmp/$sample\_aln.sai.sam
			echo "bwa samse -r "@RG\tID:$sample\tLB:$sample\tPL:ILLUMINA\tSM:$sample" $ref $tmp/$sample\_aln.sai $left > $tmp/$sample\_aln.sai.sam"

		else
			bwa aln -t 16 $ref $right > $tmp/$sample2\_aln.sai
			echo "bwa aln -t 16 $ref $right > $tmp/$sample2\_aln.sai"
			bwa sampe -r "@RG\tID:$sample\tLB:$sample\tPL:ILLUMINA\tSM:$sample" $ref $tmp/$sample\_aln.sai $tmp/$sample2\_aln.sai $left $right >  $tmp/$sample\_aln.sai.sam
			echo "bwa sampe -r "@RG\tID:$sample\tLB:$sample\tPL:ILLUMINA\tSM:$sample" $ref $tmp/$sample\_aln.sai $tmp/$sample2\_aln.sai $left $right >  $tmp/$sample\_aln.sai.sam"

                fi



		# 2. sort aln sam and convert sam to bam
		input=$tmp/$sample\_aln.sai.sam
		$picard_run/SortSam.jar I=$input O=$input.sort SO=coordinate VALIDATION_STRINGENCY=LENIENT
		echo "$picard_run/SortSam.jar I=$input O=$input.sort SO=coordinate VALIDATION_STRINGENCY=LENIENT"

		# 3. delete duplicates
		input=$input.sort
		$picard_run/MarkDuplicates.jar I=$input O=$input.dedup METRICS_FILE=$input.metrics.txt VALIDATION_STRINGENCY=LENIENT
		echo "$picard_run/MarkDuplicates.jar I=$input O=$input.dedup METRICS_FILE=$input.metrics.txt VALIDATION_STRINGENCY=LENIENT"


		# 4. add RG to distinguish sample from different source, this step can be avoid if bwa aln has added RG to results
		# addRG_bam=$dedup_bam.addRG_bam
		# $picard_run/AddOrReplaceReadGroups.jar I=$dedup_bam O=$addRG_bam SM=$sample LB=$sample ID=ILLUMINA
		input=$input.dedup
		ln -s ../$input $input.RG.bam
		echo "ln -s $input $input.RG.bam"

		# index the bam
		samtools index $input.RG.bam
		echo "samtools index $input.RG.bam"

		input=$input.RG.bam
		# 5. Indel-based realignment
		if [ "$dbsnp" == "" ]; then

			$gatk_run -nt 12 -T RealignerTargetCreator -fixMisencodedQuals -R $ref -I $input -o $input.realn.intervals # -known indels.vcf
			echo "$gatk_run -nt 12 -T RealignerTargetCreator -fixMisencodedQuals -R $ref -I $input -o $input.realn.intervals"

			#$gatk_run -nt 12 -T IndelRealigner -R $ref -I $input -targetIntervals $input.realn_interval -o $input.realn.bam

			$gatk_run -T IndelRealigner -fixMisencodedQuals -R $ref -I $input -targetIntervals $input.realn.intervals -o $input.realn.bam

			echo "$gatk_run -T IndelRealigner -R $ref -I $input -targetIntervals $input.realn.intervals -o $input.realn.bam"

			# generate the dbsnp known site vcf for non-model species
			intput=$input.realn.bam
			$gatk_run -nt 12 -T UnifiedGenotyper -fixMisencodedQuals -R $ref  -I $input -o $input.gatk_raw_vcf -stand_call_conf 30.0 -stand_emit_conf 0 -glm BOTH -rf BadCigar
			echo "$gatk_run -nt 12 -T UnifiedGenotyper -fixMisencodedQuals -R $ref  -I $input -o $input.gatk_raw_vcf -stand_call_conf 30.0 -stand_emit_conf 0 -glm BOTH -rf BadCigar"

			samtools mpileup -DSugf $ref $input | bcftools view -Ncvg - > $input.samtools_raw_vcf
			echo "samtools mpileup -DSugf $ref $input | bcftools view -Ncvg - > $input.samtools_raw_vcf"

			$gatk_run -nt 12 -T SelectVariants -R $ref --variant $input.gatk_raw_vcf --concordance $input.samtools_raw_vcf -o $input.concordance_raw_vcf
			echo "$gatk_run -nt 12 -T SelectVariants -R $ref --variant $input.gatk_raw_vcf --concordance $input.samtools_raw_vcf -o $input.concordance_raw_vcf"

			MEANQUAL=`awk '{ if ($1 !~ /#/) { total += $6; count++ } } END { print total/count }' $input.concordance_raw_vcf`

			$gatk_run -T VariantFiltration -R $ref  --filterExpression " QD < 20.0 || ReadPosRankSum < -8.0 || FS > 10.0 || QUAL < $MEANQUAL" --filterName LowQualFilter --missingValuesInExpressionsShouldEvaluateAsFailing --variant $input.concordance_raw_vcf --logging_level ERROR -o $input.concordance_flt_vcf

			echo "$gatk_run -nt 12 -T VariantFiltration -R $ref  --filterExpression " QD < 20.0 || ReadPosRankSum < -8.0 || FS > 10.0 || QUAL < $MEANQUAL" --filterName LowQualFilter --missingValuesInExpressionsShouldEvaluateAsFailing --variant $input.concordance_raw_vcf --logging_level ERROR -o $input.concordance_flt_vcf"

			# add break
			#exit 1

			grep -v Filter $input.concordance_flt_vcf > $input.known_site_vcf
			echo "grep -v Filter $input.concordance_flt_vcf > $input.known_site_vcf"
			dbsnp=$input.known_site_vcf
		fi

		realn_interval=$addRG_bam.realn.intervals
		$gatk_run -nt 12 -T RealignerTargetCreator -R $ref -I $input -o $input.realn.intervals -known $dbsnp
		echo "$gatk_run -nt 12 -T RealignerTargetCreator -R $ref -I $input -o $input.realn.intervals -known $dbsnp"


		realn_bam=$addRG_bam.realn.bam
		#$gatk_run -nt 12 -T IndelRealigner -R $ref -I $input -targetIntervals $input.realn.intervals -o $input.realn.bam -known $dbsnp

		$gatk_run -T IndelRealigner -fixMisencodedQuals -R $ref -I $input -targetIntervals $input.realn.intervals -o $input.realn.bam -known $dbsnp


		echo "$gatk_run -T IndelRealigner -fixMisencodedQuals -R $ref -I $input -targetIntervals $input.realn.intervals -o $input.realn.bam -known $dbsnp"


		# 6. Base quality score recalibration
		input=$input.realn.bam
		$gatk_run -T BaseRecalibrator -R $ref -I $input -o $input.recal_tab -knownSites $dbsnp
		echo "$gatk_run -T BaseRecalibrator -R $ref -I $input -o $input.recal_tab -knownSites $dbsnp"

		recal_bam=$realn_bam.recal.bam
		$gatk_run -T PrintReads -R $ref -I $input -BQSR $input.recal_tab -o $input.recal.bam
		echo "$gatk_run -T PrintReads -R $ref -I $input -BQSR $input.recal_tab -o $input.recal.bam"

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



