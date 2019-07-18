#!/bikn/sh
# input the list of reads

A=''
for i in `cat $1`
	do
		ARR=$i
		PE=${#ARR[@]}
		if [ "$PE" == "1" ]; then
			echo $PE
		fi
		A=`echo $A $i`

	done

echo $A
