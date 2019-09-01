#!/bin/bash


#module load parallel/20150622
flag=0
#for i in `find ./ -name \*.bam | grep \Resequencing | grep -v mapp`
#for i in `find ./ -name \*.unmapped.bam | grep \Resequencing`
for i in `find ./ -name \*.unmapped.bam.fq.fa | grep \Resequencing`
do
	file=$i\_idba_output/scaffold.fa
	#ls -sh $i.unmapped.bam
	if [ -f "$file" ]; then
    		echo "$file exist"
		continue;
	fi

	cat header.sh > run_hyalite_$flag\.sh
	echo cd $PWD >> run_hyalite_$flag\.sh

	# filter bam
	#echo "/usr/bin/time -v python ../../camelina_pangenome_analysis/sam_flt.py $i > $i\.log" >> run_hyalite_$flag\.sh
	#j=$i\.unmapped.bam

	# extract fasta from filtered bam
	j=$i
	ls -sh $j
	#echo "/usr/bin/time -v bamToFastq -i $j -fq $j\.fq &>> $i\.log" >> run_hyalite_$flag\.sh
	#echo "/usr/bin/time -v fq2fa $j\.fq $j\.fq.fa &>> $i\.log" >> run_hyalite_$flag\.sh

	# assembly
	echo "/usr/bin/time -v idba --num_threads 24 -r $j --pre_correction --mink 101 --maxk 131 -o $i\_idba_output &>> $i\_idba.log" >> run_hyalite_$flag\.sh


	# submit work
	sbatch run_hyalite_$flag\.sh
	flag=`echo $flag+1 | bc`
	echo $flag

done



