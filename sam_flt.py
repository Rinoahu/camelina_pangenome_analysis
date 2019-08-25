#!usr/bin/env python
# usage:
# python this.py foo.bam

import sys
import pysam

try:
	qry = sys.argv[1]
except:
	print('python this.py foo.bam')
	raise SystemExit()


# filter and write
f = pysam.AlignmentFile(qry, "rb")
_o1 = pysam.AlignmentFile(qry+'.mapped.bam', "wb", template=f)
_o = pysam.AlignmentFile('-', "wb", template=f)


for i in f:
	#if len(i.cigar) == 1 and i.cigar[0][0] == 0:
	#if sum([elem[1] for elem in i.cigar if elem[0]!=0]) <= 10:
	flag = sum([elem[1] for elem in i.cigar if elem[0]!=0])
	if flag < 10:
		_o1.write(i)
	else:
		_o.write(i)


f.close()
_o.close()

