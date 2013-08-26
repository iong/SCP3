function x = MatousekAffineOwen(dim, skip, Npoints)
	P = sobolset(dim, 'Skip', skip);
	P = scramble(P, 'MatousekAffineOwen');
	x = net(P, Npoints);
end
