#!/usr/bin/env python
import sys
from numpy import *
import h5py

fin=open(sys.argv[1])
fname_out = sys.argv[1][:-4]+'.h5'

row = fromstring(fin.readline(), float32, sep=' ')
N = len(row)

fout = h5py.File(fname_out, 'w')
nchunk_rows = (64*1024 / (4 * N) + 1)
dset = fout.create_dataset('M', (N,N), '<f4', chunks=(nchunk_rows,N), compression='gzip')

dset[0,:] = row
for i in range(1,N):
    dset[i,:] = fromstring(fin.readline(), float32, sep=' ')

fin.close()
fout.close()
