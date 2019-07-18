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
#$picard_run/CreateSequenceDictionary.jar R=$ref O=$ref.dict
#echo "$picard_run/CreateSequenceDictionary.jar R=$ref O=$ref.dict"
echo "gatk CreateSequenceDictionary R=$ref O=$ref.dict"
gatk CreateSequenceDictionary R=$ref O=$ref.dict


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
		if [ "$right" == "" ]; then
			echo "bwa mem -t 16 -R "@RG\tID:$sample\tLB:$sample\tPL:ILLUMINA\tSM:$sample" $ref $left > $tmp/$sample\_aln.sai.sam"
			bwa mem -t 16 -R "@RG\tID:$sample\tLB:$sample\tPL:ILLUMINA\tSM:$sample" $ref $left > $tmp/$sample\_aln.sai.sam

		else
			echo "bwa mem -t 16 -R "@RG\tID:$sample\tLB:$sample\tPL:ILLUMINA\tSM:$sample" $ref $left $right >  $tmp/$sample\_aln.sai.sam"
			bwa mem -t 16 -R "@RG\tID:$sample\tLB:$sample\tPL:ILLUMINA\tSM:$sample" $ref $left $right >  $tmp/$sample\_aln.sai.sam
        fi


		# 2. sort aln sam and convert sam to bam
		input=$tmp/$sample\_aln.sai.sam
		#$picard_run/SortSam.jar I=$input O=$input.sort SO=coordinate VALIDATION_STRINGENCY=LENIENT
		#echo "$picard_run/SortSam.jar I=$input O=$input.sort SO=coordinate VALIDATION_STRINGENCY=LENIENT"
		echo "gatk SortSam I=$input O=$input.sort SO=coordinate VALIDATION_STRINGENCY=LENIENT"
		gatk SortSam I=$input O=$input.sort SO=coordinate VALIDATION_STRINGENCY=LENIENT


		# 3. delete duplicates
		input=$input.sort
		#$picard_run/MarkDuplicates.jar I=$input O=$input.dedup METRICS_FILE=$input.metrics.txt VALIDATION_STRINGENCY=LENIENT
		#echo "$picard_run/MarkDuplicates.jar I=$input O=$input.dedup METRICS_FILE=$input.metrics.txt VALIDATION_STRINGENCY=LENIENT"

        echo "gatk MarkDuplicates I=$input O=$input.dedup METRICS_FILE=$input.metrics.txt VALIDATION_STRINGENCY=LENIENT"
        gatk MarkDuplicates I=$input O=$input.dedup METRICS_FILE=$input.metrics.txt VALIDATION_STRINGENCY=LENIENT


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

			echo "gatk HaplotypeCaller -R $ref  -I $input -o $input.gatk_raw_vcf -stand_call_conf 30.0 --emit-ref-confidence GVCF"
            gatk HaplotypeCaller -R $ref  -I $input -o $input.gatk_raw_vcf -stand_call_conf 30.0 --emit-ref-confidence GVCF

			echo "samtools mpileup -DSugf $ref $input | bcftools view -Ncvg - > $input.samtools_raw_vcf"
			samtools mpileup -DSugf $ref $input | bcftools view -Ncvg - > $input.samtools_raw_vcf

			#$gatk_run -nt 12 -T SelectVariants -R $ref --variant $input.gatk_raw_vcf --concordance $input.samtools_raw_vcf -o $input.concordance_raw_vcf
			echo "gatk SelectVariants -select-type SNP -R $ref --variant $input.gatk_raw_vcf --concordance $input.samtools_raw_vcf -o $input.concordance_raw_vcf"
			gatk SelectVariants -select-type SNP -R $ref --variant $input.gatk_raw_vcf --concordance $input.samtools_raw_vcf -o $input.concordance_raw_snp_vcf

			echo "gatk SelectVariants -select-type INDEL -R $ref --variant $input.gatk_raw_vcf --concordance $input.samtools_raw_vcf -o $input.concordance_raw_vcf"
			gatk SelectVariants -select-type INDEL -R $ref --variant $input.gatk_raw_vcf --concordance $input.samtools_raw_vcf -o $input.concordance_raw_indel_vcf

            echo "gatk MergeVcfs -I $input.concordance_raw_snp_vcf -I $input.concordance_raw_indel_vcf -O $input.concordance_raw_vcf"
            gatk MergeVcfs -I $input.concordance_raw_snp_vcf -I $input.concordance_raw_indel_vcf -O $input.concordance_raw_vcf

			MEANQUAL=`awk '{ if ($1 !~ /#/) { total += $6; count++ } } END { print total/count }' $input.concordance_raw_vcf`

			echo "gatk VariantFiltration --filter-expression " QD < 20.0 || ReadPosRankSum < -8.0 || FS > 10.0 || QUAL < $MEANQUAL" --filterName LowQualFilter --missing-values-evaluate-as-failing true --variant $input.concordance_raw_vcf --verbosity ERROR -R $ref -o $input.concordance_flt_vcf"

            gatk VariantFiltration --filter-expression " QD < 20.0 || ReadPosRankSum < -8.0 || FS > 10.0 || QUAL < $MEANQUAL" --filterName LowQualFilter --missing-values-evaluate-as-failing true --variant $input.concordance_raw_vcf --verbosity ERROR -R $ref -o $input.concordance_flt_vcf

			# add break
			#exit 1

			grep -v Filter $input.concordance_flt_vcf > $input.known_site_vcf
			echo "grep -v Filter $input.concordance_flt_vcf > $input.known_site_vcf"
			dbsnp=$input.known_site_vcf
		fi



		# 6. Base quality score recalibration
        echo "gatk BaseRecalibrator -R $ref -I $input -o $input.recal_tab -knownSites $dbsnp"
		gatk BaseRecalibrator -R $ref -I $input -O $input.recal_tab --known-sites $dbsnp
        gatk ApplyBQSR -R $ref -I $input --bqsr-recal-file $input.recal_tab -O $input.recal.bam



        exit 1
		recal_bam=$realn_bam.recal.bam
        echo "gatk PrintReads -R $ref -I $input -BQSR $input.recal_tab -o $input.recal.bam"
		gatk PrintReads -R $ref -I $input -BQSR $input.recal_tab -o $input.recal.bam


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



