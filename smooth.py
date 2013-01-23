from numpy import *
def smooth(omega, mu, w):
	grid=arange(8000)
	s = grid*0.0
	for (omega0, mu0) in zip(omega, mu):
		s+= mu0*exp(-(grid - omega0)**2/(2*w**2))
	s/= w*sqrt(2.0*pi)
	return (grid, s)
