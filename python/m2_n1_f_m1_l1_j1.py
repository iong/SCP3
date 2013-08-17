def m2_n1_f_m1_l1_j1(q, V, Vq, m, n, l, j):
	return (4.0*V*q[:,j]*q[:,l]*q[:,m]*q[:,n] + (q[:,l]*q[:,m]*q[:,n]*Vq[:,j] + q[:,j]*(q[:,m]*q[:,n]*Vq[:,l] + q[:,l]*q[:,n]*Vq[:,m] + q[:,l]*q[:,m]*Vq[:,n]))*(-1.0 + 2.0*sqr(q[:,m])))/(2.*sqrt(2.0))
