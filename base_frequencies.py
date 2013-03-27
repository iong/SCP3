def basis_frequencies(w0, Nmodes2, Nmodes3):

	Nmodes = len(w0)

	Nstates = 1 + Nmodes + Nmodes2*(Nmodes2 + 1)/2
        Nstates += (Nmodes3 * (1 + Nmodes3) * (2 + Nmodes3))/6

	wout = zeros((Nstates,))
	C = 1

	wout[C : C + Nmodes] = w0
	C += Nmodes

	wout[C : C + Nmodes2] = 2*w0[:Nmodes2]
	C += Nmodes2

	for k in range(Nmodes2):
		M = Nmodes2 - k - 1
		wout[C : C + M] =  w0[k+1 : Nmodes2] + w0[k]
		C += M
	
	wout[C : C + Nmodes3] = 3*w0[:Nmodes3]
	C += Nmodes3

	for k in range(Nmodes3):
		M = Nmodes3 - k - 1
		wout[C   : C + 2*M : 2] =   w0[k+1 : Nmodes3] + 2*w0[k]
		wout[C+1 : C + 2*M : 2] = 2*w0[k+1 : Nmodes3] +   w0[k]
		C += 2*M

	for k in range(Nmodes3):
		for l in range(k+1, Nmodes3):
			M = Nmodes3 - l - 1
			wout[C : C + M] = w0[l+1:Nmodes3] + w0[k] + w0[l]
			C += M
	
	print C, len(wout)

	return wout


def state_list(Nmodes, Nmodes2, Nmodes3):

	sl = [(-1,)]
	sl += [(x,) for x in range(Nmodes)] + [(x,x) for x in range(Nmodes2)]

	for k in range(Nmodes2):
		sl += [(k, l) for l in range(k+1,Nmodes2)]
	
	sl += [(x, x, x) for x in range(Nmodes3)]

	for k in range(Nmodes3):
		for l in range(k+1,Nmodes3):
			sl += [(k, k, l), (k, l, l)]

	for k in range(Nmodes3):
		for l in range(k+1, Nmodes3):
			sl += [(k, l, j) for j in range(l+1,Nmodes3)]
	
	return sl
