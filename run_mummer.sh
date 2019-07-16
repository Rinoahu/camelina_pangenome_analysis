#!/bin/sh
if [ $# -lt 2 ]
then
	echo 'usage:'
	echo 'this_script.sh ref.fsa qry.fsa'
	exit 1

fi


ref=$1
qry=$2
prog=$3
echo 'query is', $qry
echo 'ref is', $ref

if [[ "$prog" = "" ]]
	then prog='promer'
fi

echo 'prog is', $prog


#promer --prefix=ref_qry $ref $qry --mum
#$prog --prefix=ref_qry $ref $qry --mum


#show-coords -rcl ref_qry.delta > ref_qry.coords

#show-aligns ref_qry.delta $ref $qry > ref_qry.aligns

delta-filter -q -r ref_qry.delta -i 98 -l 5000 > ref_qry.filter
#delta-filter -q -r ref_qry.delta -i 10 -l 30 > ref_qry.filter

mummerplot ref_qry.filter -R $ref -Q $qry --layout --postscript --prefix=ref_qry.filter
