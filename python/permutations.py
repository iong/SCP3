def perm(f, fcut):
	n = floor(fcut/f)

	f = flipud(sort(f))

	budget = 0.0*f

	for i0 in range(int(fcut/f[0]+1.0)):
		budget[0] = fcut - i0*f[1]
		for i1 in range(int(budget[0]/f[1] + 1.0)):
			budget[1] = budget[0] - i1*f[1]
			for i2 in range(int(budget[1]/f[2]  +1.0)):
				budget[2] = budget[1] - i2*f[2]
				print i2, i1, i0

def rperm_(f, n, budget, level, lNperms):
	for n[level] in range(int(budget[level-1]/f[level]+1.0)):
		budget[level] = budget[level - 1] - n[level]*f[level]
		if level+1 < len(f):
			rperm_(f, n, budget, level+1, lNperms)
		else:
			lNperms += [budget[-1] - budget[level]]

def rperm(f, fcut):
	n=zeros(f.shape)
	budget=append(0.0*f, fcut)
	lNperms=[0]
	rperm_(f, n, budget, 0, lNperms)
	return lNperms
