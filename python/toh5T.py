#!/usr/bin/env python
import sys
from numpy import *
import h5py

fname_in=sys.argv[1]
fname_out = sys.argv[1][:-4]+'.h5'

M = loadtxt(fname_in, dtype=float32)
N = M.shape[0]

fout = h5py.File(fname_out, 'w')
nchunk_rows = (64*1024 / (4 * N) + 1)
dset = fout.create_dataset('M', (N,N), '<f4', chunks=(nchunk_rows, N), compression='gzip')

dset[:,:] = M.T
fout.close()
