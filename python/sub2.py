#!/usr/bin/python

# Incomplete script for processing TripleExcitations.h
# Kept for reference.

f.seek(0)
fout=open('TE.h', 'w')

a=re.compile("square\((\w+)\[i\]\)")

for l in f:
    vars = set(a.findall(l))
    if vars:
        print vars
        for v in vars:
            fout.write("        double _sqr_%si = %s[i]*%s[i];\n"%(v,v,v))
        fout.write(a.sub("_sqr_\\1i", l))
    else:
        fout.write(l)

fout.close()

