#!/usr/bin/python

# Incomplete script for processing TripleExcitations.h
# Kept for reference.

import re

f.seek(0)
fout=open('TE.h', 'w')

for l in f:
    fout.write(l)
    if a.match(l):
        fout.write("{\n")
        vars=b.findall(l)
        for v in vars:
            fout.write("    double *Vq%s = __assume_aligned(Vq.colptr(%s), 32);\n"%(v,v))
            fout.write("    double * q%s = __assume_aligned( q.colptr(%s), 32);\n"%(v,v))
        l = f.next()

fout.close()

