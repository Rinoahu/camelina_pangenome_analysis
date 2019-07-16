#!usr/bin/env python

import sys
from Bio import SeqIO
import re

p = re.compile(r'[N|n]+')

try:
    qry = sys.argv[1]
except:
    print('python this.py foo.fsa')
    raise SystemExit()

seqs = SeqIO.parse(qry, 'fasta')

lns = []
Ns = 0
for i in seqs:
    seq = str(i.seq)
    Ns += sum(map(len, p.findall(seq)))
    lns.append(len(seq))

lns.sort(reverse=True)
# find N50
N50 = sum(lns) // 2
n50 = 0
idx = 0
for i in xrange(len(lns)):
    if n50 >= N50:
        idx = i
        #print('idx is', idx, n50)
        break
    else:
        n50 += lns[i]

print('genome size', sum(lns))
print('   # of scf', len(lns))
print('    max scf', lns[0])
print('    min scf', lns[-1])
print('     median', lns[len(lns)//2])
print('        N50', lns[idx])
print('         Ns', Ns)
