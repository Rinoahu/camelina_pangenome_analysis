#!/bin/bash

j=""
for i in merge_*.vcf
do
	j="$j $i"
	#echo $j
	bcftools index $i
done

bcftools concat merge_0.vcf merge_10.vcf merge_11.vcf merge_12.vcf merge_13.vcf merge_14.vcf merge_15.vcf merge_16.vcf merge_17.vcf merge_18.vcf merge_19.vcf merge_1.vcf merge_20.vcf merge_21.vcf merge_22.vcf merge_23.vcf merge_24.vcf merge_25.vcf merge_26.vcf merge_27.vcf merge_28.vcf merge_29.vcf merge_2.vcf merge_30.vcf merge_31.vcf merge_32.vcf merge_33.vcf merge_34.vcf merge_35.vcf merge_36.vcf merge_37.vcf merge_38.vcf merge_39.vcf merge_3.vcf merge_40.vcf merge_41.vcf merge_42.vcf merge_4.vcf merge_5.vcf merge_6.vcf merge_7.vcf merge_81.vcf merge_82.vcf merge_83.vcf merge_84.vcf merge_85.vcf merge_86.vcf merge_87.vcf merge_88.vcf merge_89.vcf merge_8.vcf merge_90.vcf merge_91.vcf merge_92.vcf merge_93.vcf merge_94.vcf merge_9.vcf -o merge_all.bcf -a -D
