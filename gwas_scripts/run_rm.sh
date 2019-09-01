#!/bin/bash


#module load parallel/20150622
flag=0
#for i in `find ./ -name \*.bam | grep \Resequencing | grep -v mapp`
#for i in `find ./ -name \*.unmapped.bam | grep \Resequencing`
for i in `find ./ -name \*.unmapped.bam.fq | grep \Resequencing`
do
	file=$i\.fa\_idba_output/scaffold.fa
	#ls -sh $i.unmapped.bam
	if [ -f "$file" ]; then
    		echo "$file exist"
		rm $i
		rm $i\.fa
		#rm $i\_idba_output/align-*
		#rm $i\_idba_output/contig-*.fa
		#rm $i\_idba_output/graph-*.fa
		#rm $i\_idba_output/kmer	
		continue;
	fi


done



