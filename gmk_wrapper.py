#!usr/bin/env python
from Bio import SeqIO
import os, sys
try:
	qry, mod = sys.argv[1:3]
except:
	print('python this.py foo.fsa foo.mod')
	raise SystemExit()

seqs = SeqIO.parse(qry, 'fasta')

for seq in seqs:
	_o = open(qry + '_tmp.fsa', 'w')
	SeqIO.write([seq], _o, 'fasta')
	_o.close()

	cmd = 'gmhmme3 -m %s -f gtf -o %s_tmp.fsa.gff3 %s_tmp.fsa'
	os.system(cmd%(mod, qry, qry))

	if not os.path.isfile('%s_tmp.fsa.gff3'%qry):
		continue

	qid = seq.id
	#print('qid', qid)
	f = open('%s_tmp.fsa.gff3'%qry, 'r')
	for i in f:
		j = i[:-1].split('\t')
		j[0] = j[0] == 'seq' and qid or j[0]
		print('\t'.join(j))

	f.close()
	os.system('rm %s_tmp.fsa.gff3'%qry)
