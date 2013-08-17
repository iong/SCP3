from numpy import *
def smooth(omega, mu, w):
	grid = arange(8000)
	s = grid*0.0
#	do i=1,size(omega)
#		omega0=omega(i)
#		mu0=mu(i)
#		do j=1,8000
#			s(j) += mu0*exp(-(grid(j) - omega0)**2/(2*w**2))
	for (omega0, mu0) in zip(omega, mu):
		s += mu0*exp(-(grid - omega0)**2/(2*w**2))
	s /= w*sqrt(2.0*pi)
	return (grid, s)
