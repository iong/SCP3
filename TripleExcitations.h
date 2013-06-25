double f_k1(int k)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += Vq(i, k);
	}
	return ret / sqrt(2.0);
}
double f_k2(int k)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += q(i, k)*Vq(i, k);
	}
	return ret / sqrt(2.);
}
double f_k3(int k)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += 4.*V(i)*q(i, k) + Vq(i, k)*(-3. + 2.*square(q(i, k)));
	}
	return ret / (2.*sqrt(3.));
}
double f_k1_l1(int k, int l)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += q(i, l)*Vq(i, k) + q(i, k)*Vq(i, l);
	}
	return ret / 2.;
}
double f_k1_l2(int k, int l)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += q(i, k)*q(i, l)*Vq(i, l) + Vq(i, k)*(-0.5 + square(q(i, l)));
	}
	return ret / 2.;
}
double f_k2_l1(int k, int l)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += q(i, k)*q(i, l)*Vq(i, k) + Vq(i, l)*(-0.5 + square(q(i, k)));
	}
	return ret / 2.;
}
double f_k1_l1_j1(int k, int l, int j)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += sqrt(2.)*(q(i, k)*q(i, l)*Vq(i, j) + q(i, j)*q(i, l)*Vq(i, k) + q(i, j)*q(i, k)*Vq(i, l));
	}
	return ret / 3.;
}

double m1_f_k1(int m, int k)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += q(i, m)*Vq(i, k) + q(i, k)*Vq(i, m);
	}
	return ret / 2.;
}

double m1_f_k2(int m, int k)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += q(i, k)*q(i, m)*Vq(i, k) + Vq(i, m)*(-0.5 + square(q(i, k)));
	}
	return ret / 2.;
}


double m1_f_k3(int m, int k)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += 4.*V(i)*q(i, k)*q(i, m) + (q(i, m)*Vq(i, k) + q(i, k)*Vq(i, m))*(-3. + 2.*square(q(i, k)));
	}
	return ret / (2.*sqrt(6.));
}


double m1_f_k1_l1(int m, int k, int l)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += sqrt(2.)*(q(i, l)*q(i, m)*Vq(i, k) + q(i, k)*q(i, m)*Vq(i, l) + q(i, k)*q(i, l)*Vq(i, m));
	}
	return ret / 3.;
}


double m1_f_k1_l2(int m, int k, int l)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += 2.*q(i, k)*q(i, l)*q(i, m)*Vq(i, l) + q(i, m)*Vq(i, k)*(-1. + 2.*square(q(i, l))) + q(i, k)*Vq(i, m)*(-1. + 2.*square(q(i, l)));
	}
	return ret / (3.*sqrt(2.));
}
double m1_f_k2_l1(int m, int k, int l)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += 2.*q(i, k)*q(i, l)*q(i, m)*Vq(i, k) + (q(i, m)*Vq(i, l) + q(i, l)*Vq(i, m))*(-1. + 2.*square(q(i, k)));
	}
	return ret / (3.*sqrt(2.));
}
double m1_f_k1_l1_j1(int m, int k, int l, int j)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += q(i, k)*q(i, l)*q(i, m)*Vq(i, j) + q(i, j)*(q(i, l)*q(i, m)*Vq(i, k) + q(i, k)*q(i, m)*Vq(i, l) + q(i, k)*q(i, l)*Vq(i, m));
	}
	return ret / 2.;
}
double m2_f_k2(int m, int k)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += q(i, m)*Vq(i, m)*(-1. + 2.*square(q(i, k))) + q(i, k)*Vq(i, k)*(-1. + 2.*square(q(i, m)));
	}
	return ret / 4.;
}
double m2_f_k3(int m, int k)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += 4.*V(i)*q(i, k)*(-1. + 2.*square(q(i, m))) + (-3. + 2.*square(q(i, k)))*(2.*q(i, k)*q(i, m)*Vq(i, m) + Vq(i, k)*(-1. + 2.*square(q(i, m))));
	}
	return ret / (4.*sqrt(6.));
}
double m2_f_k1_l1(int m, int k, int l)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += 2.*q(i, k)*q(i, l)*q(i, m)*Vq(i, m) + (q(i, l)*Vq(i, k) + q(i, k)*Vq(i, l))*(-1. + 2.*square(q(i, m)));
	}
	return ret / (3.*sqrt(2.));
}
double m2_f_k1_l2(int m, int k, int l)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += 2.*q(i, k)*q(i, m)*Vq(i, m)*(-1. + 2.*square(q(i, l))) + (2.*q(i, k)*q(i, l)*Vq(i, l) + Vq(i, k)*(-1. + 2.*square(q(i, l))))*(-1. + 2.*square(q(i, m)));
	}
	return ret / (6.*sqrt(2.));
}
double m2_f_k2_l1(int m, int k, int l)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += 2.*q(i, k)*q(i, l)*Vq(i, k)*(-1. + 2.*square(q(i, m))) + (-1. + 2.*square(q(i, k)))*(2.*q(i, l)*q(i, m)*Vq(i, m) + Vq(i, l)*(-1. + 2.*square(q(i, m))));
	}
	return ret / (6.*sqrt(2.));
}
double m2_f_k1_l1_j1(int m, int k, int l, int j)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += 2.*q(i, j)*q(i, k)*q(i, l)*q(i, m)*Vq(i, m) + (q(i, k)*q(i, l)*Vq(i, j) + q(i, j)*q(i, l)*Vq(i, k) + q(i, j)*q(i, k)*Vq(i, l))*(-1. + 2.*square(q(i, m)));
	}
	return ret / 4.;
}
double m3_f_k3(int m, int k)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += 4.*V(i)*q(i, k)*q(i, m)*(-3. + 2.*square(q(i, m))) + (-3. + 2.*square(q(i, k)))*(4.*V(i)*q(i, k)*q(i, m) + (q(i, m)*Vq(i, k) + q(i, k)*Vq(i, m))*(-3. + 2.*square(q(i, m))));
	}
	return ret / 12.;
}
double m3_f_k1_l1(int m, int k, int l)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += 4.*V(i)*q(i, k)*q(i, l)*q(i, m) + (q(i, l)*q(i, m)*Vq(i, k) + q(i, k)*q(i, m)*Vq(i, l) + q(i, k)*q(i, l)*Vq(i, m))*(-3. + 2.*square(q(i, m)));
	}
	return ret / (3.*sqrt(3.));
}
double m3_f_k1_l2(int m, int k, int l)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += 4.*V(i)*q(i, k)*q(i, m)*(-1. + 2.*square(q(i, l))) + (2.*q(i, k)*q(i, l)*q(i, m)*Vq(i, l) + q(i, m)*Vq(i, k)*(-1. + 2.*square(q(i, l))) + q(i, k)*Vq(i, m)*(-1. + 2.*square(q(i, l))))*(-3. + 2.*square(q(i, m)));
	}
	return ret / (6.*sqrt(3.));
}
double m3_f_k2_l1(int m, int k, int l)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += 4.*V(i)*q(i, l)*q(i, m)*(-1. + 2.*square(q(i, k))) + (2.*q(i, k)*q(i, l)*q(i, m)*Vq(i, k) + (q(i, m)*Vq(i, l) + q(i, l)*Vq(i, m))*(-1. + 2.*square(q(i, k))))*(-3. + 2.*square(q(i, m)));
	}
	return ret / (6.*sqrt(3.));
}
double m3_f_k1_l1_j1(int m, int k, int l, int j)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += 4.*V(i)*q(i, j)*q(i, k)*q(i, l)*q(i, m) + (q(i, k)*q(i, l)*q(i, m)*Vq(i, j) + q(i, j)*(q(i, l)*q(i, m)*Vq(i, k) + q(i, k)*q(i, m)*Vq(i, l) + q(i, k)*q(i, l)*Vq(i, m)))*(-3. + 2.*square(q(i, m)));
	}
	return ret / (2.*sqrt(6.));
}
double m1_n1_f_k1_l1(int m, int n, int k, int l)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += q(i, l)*q(i, m)*q(i, n)*Vq(i, k) + q(i, k)*(q(i, m)*q(i, n)*Vq(i, l) + q(i, l)*q(i, n)*Vq(i, m) + q(i, l)*q(i, m)*Vq(i, n));
	}
	return ret / 2.;
}
double m1_n1_f_k1_l2(int m, int n, int k, int l)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += 2.*q(i, k)*q(i, l)*q(i, m)*q(i, n)*Vq(i, l) + q(i, m)*q(i, n)*Vq(i, k)*(-1. + 2.*square(q(i, l))) + q(i, k)*(q(i, n)*Vq(i, m) + q(i, m)*Vq(i, n))*(-1. + 2.*square(q(i, l)));
	}
	return ret / 4.;
}
double m1_n1_f_k2_l1(int m, int n, int k, int l)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += 2.*q(i, k)*q(i, l)*q(i, m)*q(i, n)*Vq(i, k) + (q(i, m)*q(i, n)*Vq(i, l) + q(i, l)*q(i, n)*Vq(i, m) + q(i, l)*q(i, m)*Vq(i, n))*(-1. + 2.*square(q(i, k)));
	}
	return ret / 4.;
}
double m1_n1_f_k1_l1_j1(int m, int n, int k, int l, int j)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += 2.*sqrt(2.)*(q(i, k)*q(i, l)*q(i, m)*q(i, n)*Vq(i, j) + q(i, j)*(q(i, l)*q(i, m)*q(i, n)*Vq(i, k) + q(i, k)*(q(i, m)*q(i, n)*Vq(i, l) + q(i, l)*q(i, n)*Vq(i, m) + q(i, l)*q(i, m)*Vq(i, n))));
	}
	return ret / 5.;
}
double m1_n2_f_k1_l2(int m, int n, int k, int l)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += 2.*q(i, k)*q(i, m)*q(i, n)*Vq(i, n)*(-1. + 2.*square(q(i, l))) + (2.*q(i, k)*q(i, l)*q(i, m)*Vq(i, l) + q(i, m)*Vq(i, k)*(-1. + 2.*square(q(i, l))) + q(i, k)*Vq(i, m)*(-1. + 2.*square(q(i, l))))*(-1. + 2.*square(q(i, n)));
	}
	return ret / 8.;
}
double m1_n2_f_k2_l1(int m, int n, int k, int l)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += 2.*q(i, k)*q(i, l)*q(i, m)*Vq(i, k)*(-1. + 2.*square(q(i, n))) + (-1. + 2.*square(q(i, k)))*(2.*q(i, l)*q(i, m)*q(i, n)*Vq(i, n) + (q(i, m)*Vq(i, l) + q(i, l)*Vq(i, m))*(-1. + 2.*square(q(i, n))));
	}
	return ret / 8.;
}
double m1_n2_f_k1_l1_j1(int m, int n, int k, int l, int j)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += sqrt(2.)*(2.*q(i, j)*q(i, k)*q(i, l)*q(i, m)*q(i, n)*Vq(i, n) + (q(i, k)*q(i, l)*q(i, m)*Vq(i, j) + q(i, j)*(q(i, l)*q(i, m)*Vq(i, k) + q(i, k)*q(i, m)*Vq(i, l) + q(i, k)*q(i, l)*Vq(i, m)))*(-1. + 2.*square(q(i, n))));
	}
	return ret / 5.;
}
double m2_n1_f_k2_l1(int m, int n, int k, int l)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += 2.*q(i, k)*q(i, l)*q(i, n)*Vq(i, k)*(-1. + 2.*square(q(i, m))) + (-1. + 2.*square(q(i, k)))*(2.*q(i, l)*q(i, m)*q(i, n)*Vq(i, m) + q(i, n)*Vq(i, l)*(-1. + 2.*square(q(i, m))) + q(i, l)*Vq(i, n)*(-1. + 2.*square(q(i, m))));
	}
	return ret / 8.;
}
double m2_n1_f_k1_l1_j1(int m, int n, int k, int l, int j)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += sqrt(2.)*(2.*q(i, j)*q(i, k)*q(i, l)*q(i, m)*q(i, n)*Vq(i, m) + q(i, n)*(q(i, k)*q(i, l)*Vq(i, j) + q(i, j)*q(i, l)*Vq(i, k) + q(i, j)*q(i, k)*Vq(i, l))*(-1. + 2.*square(q(i, m))) + q(i, j)*q(i, k)*q(i, l)*Vq(i, n)*(-1. + 2.*square(q(i, m))));
	}
	return ret / 5.;
}
double m1_n1_p1_f_k1_l1_j1(int m, int n, int p, int k, int l, int j)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += 2.*(q(i, k)*q(i, l)*q(i, m)*q(i, n)*q(i, p)*Vq(i, j) + q(i, j)*q(i, l)*q(i, m)*q(i, n)*q(i, p)*Vq(i, k) + q(i, j)*q(i, k)*(q(i, m)*q(i, n)*q(i, p)*Vq(i, l) + q(i, l)*q(i, n)*q(i, p)*Vq(i, m) + q(i, l)*q(i, m)*q(i, p)*Vq(i, n) + q(i, l)*q(i, m)*q(i, n)*Vq(i, p)));
	}
	return ret / 3.;
}
double m1_f_m1(int m)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += V(i) + q(i, m)*Vq(i, m);
	}
	return ret;
}
double m1_f_m2(int m)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += 2.*V(i)*q(i, m) + Vq(i, m)*(-0.5 + square(q(i, m)));
	}
	return ret;
}
double m1_f_m3(int m)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += q(i, m)*Vq(i, m)*(-3. + 2.*square(q(i, m))) + V(i)*(-3. + 6.*square(q(i, m)));
	}
	return ret / sqrt(6.);
}
double m1_f_m1_l1(int m, int l)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += V(i)*q(i, l) + q(i, m)*(q(i, m)*Vq(i, l) + q(i, l)*Vq(i, m));
	}
	return ret / sqrt(2.);
}
double m1_f_m1_l2(int m, int l)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += (V(i) + q(i, m)*Vq(i, m))*(-1. + 2.*square(q(i, l))) + 2.*q(i, l)*Vq(i, l)*square(q(i, m));
	}
	return ret / (2.*sqrt(2.));
}
double m1_f_m2_l1(int m, int l)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += 4.*V(i)*q(i, l)*q(i, m) + (q(i, m)*Vq(i, l) + q(i, l)*Vq(i, m))*(-1. + 2.*square(q(i, m)));
	}
	return ret / (2.*sqrt(2.));
}
double m1_f_m1_l1_j1(int m, int l, int j)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += 2.*(V(i)*q(i, j)*q(i, l) + q(i, m)*(q(i, l)*q(i, m)*Vq(i, j) + q(i, j)*q(i, m)*Vq(i, l) + q(i, j)*q(i, l)*Vq(i, m)));
	}
	return ret / 3.;
}
double m2_f_m2(int m)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += q(i, m)*Vq(i, m)*(-1. + square(q(i, m))) + V(i)*(-0.5 + 3.*square(q(i, m)));
	}
	return ret;
}

double m2_f_m3(int m)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += Vq(i, m)*(3. + 4.*square(square(q(i, m))) - 8.*square(q(i, m))) + 16.*V(i)*q(i, m)*(-1. + square(q(i, m)));
	}
	return ret / (2.*sqrt(6.));
}
double m2_f_m1_l1(int m, int l)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += 4.*V(i)*q(i, l)*q(i, m) + (q(i, m)*Vq(i, l) + q(i, l)*Vq(i, m))*(-1. + 2.*square(q(i, m)));
	}
	return ret / (2.*sqrt(2.));
}
double m2_f_m1_l2(int m, int l)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += 4.*V(i)*q(i, m)*(-1. + 2.*square(q(i, l))) + (2.*q(i, l)*q(i, m)*Vq(i, l) + Vq(i, m)*(-1. + 2.*square(q(i, l))))*(-1. + 2.*square(q(i, m)));
	}
	return ret / (4.*sqrt(2.));
}
double m2_f_m2_l1(int m, int l)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += Vq(i, l) + 2.*(2.*q(i, m)*(q(i, m)*Vq(i, l) + q(i, l)*Vq(i, m))*(-1. + square(q(i, m))) + V(i)*q(i, l)*(-1. + 6.*square(q(i, m))));
	}
	return ret / (4.*sqrt(2.));
}
double m2_f_m1_l1_j1(int m, int l, int j)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += 4.*V(i)*q(i, j)*q(i, l)*q(i, m) + (q(i, l)*q(i, m)*Vq(i, j) + q(i, j)*q(i, m)*Vq(i, l) + q(i, j)*q(i, l)*Vq(i, m))*(-1. + 2.*square(q(i, m)));
	}
	return ret / 3.;
}
double m3_f_m3(int m)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += (-3. + 2.*square(q(i, m)))*(q(i, m)*Vq(i, m)*(-3. + 2.*square(q(i, m))) + V(i)*(-3. + 1.0*square(q(i, m))));
	}
	return ret / 6.;
}
double m3_f_m1_l1(int m, int l)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += q(i, m)*(q(i, m)*Vq(i, l) + q(i, l)*Vq(i, m))*(-3. + 2.*square(q(i, m))) + 3.*V(i)*q(i, l)*(-1. + 2.*square(q(i, m)));
	}
	return ret / (2.*sqrt(3.));
}
double m3_f_m1_l2(int m, int l)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += q(i, m)*(2.*q(i, l)*q(i, m)*Vq(i, l) + Vq(i, m)*(-1. + 2.*square(q(i, l))))*(-3. + 2.*square(q(i, m))) + 3.*V(i)*(-1. + 2.*square(q(i, l)))*(-1. + 2.*square(q(i, m)));
	}
	return ret / (4.*sqrt(3.));
}
double m3_f_m2_l1(int m, int l)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += (q(i, m)*Vq(i, l) + q(i, l)*Vq(i, m))*(3. + 4.*square(square(q(i, m))) - 8.*square(q(i, m))) + 16.*V(i)*q(i, l)*q(i, m)*(-1. + square(q(i, m)));
	}
	return ret / (4.*sqrt(3.));
}
double m3_f_m1_l1_j1(int m, int l, int j)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += sqrt(0.6666666666666666)*(q(i, m)*(q(i, l)*q(i, m)*Vq(i, j) + q(i, j)*q(i, m)*Vq(i, l) + q(i, j)*q(i, l)*Vq(i, m))*(-3. + 2.*square(q(i, m))) + 3.*V(i)*q(i, j)*q(i, l)*(-1. + 2.*square(q(i, m))));
	}
	return ret / 3.;
}
double m1_n1_f_m1_l1(int m, int n, int l)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += 2.*(V(i)*q(i, l)*q(i, n) + q(i, m)*(q(i, m)*q(i, n)*Vq(i, l) + q(i, l)*q(i, n)*Vq(i, m) + q(i, l)*q(i, m)*Vq(i, n)));
	}
	return ret / 3.;
}
double m1_n1_f_m1_l2(int m, int n, int l)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += V(i)*q(i, n)*(-1. + 2.*square(q(i, l))) + q(i, m)*(2.*q(i, l)*q(i, m)*q(i, n)*Vq(i, l) + (q(i, n)*Vq(i, m) + q(i, m)*Vq(i, n))*(-1. + 2.*square(q(i, l))));
	}
	return ret / 3.;
}
double m1_n1_f_m2_l1(int m, int n, int l)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += 4.*V(i)*q(i, l)*q(i, m)*q(i, n) + (q(i, m)*q(i, n)*Vq(i, l) + q(i, l)*q(i, n)*Vq(i, m) + q(i, l)*q(i, m)*Vq(i, n))*(-1. + 2.*square(q(i, m)));
	}
	return ret / 3.;
}
double m1_n1_f_m1_l1_j1(int m, int n, int l, int j)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += V(i)*q(i, j)*q(i, l)*q(i, n) + q(i, m)*(q(i, l)*q(i, m)*q(i, n)*Vq(i, j) + q(i, j)*(q(i, m)*q(i, n)*Vq(i, l) + q(i, l)*q(i, n)*Vq(i, m) + q(i, l)*q(i, m)*Vq(i, n)));
	}
	return ret / sqrt(2.);
}
double m1_n2_f_m1_l2(int m, int n, int l)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += q(i, n)*Vq(i, n)*(-1. + 2.*square(q(i, l)))*square(q(i, m)) + ((V(i) + q(i, m)*Vq(i, m))*(-1. + 2.*square(q(i, l)))*(-1. + 2.*square(q(i, n))))/2. + q(i, l)*Vq(i, l)*square(q(i, m))*(-1. + 2.*square(q(i, n)));
	}
	return ret / 3.;
}
double m1_n2_f_m2_l1(int m, int n, int l)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += 4.*V(i)*q(i, l)*q(i, m)*(-1. + 2.*square(q(i, n))) + (-1. + 2.*square(q(i, m)))*(2.*q(i, l)*q(i, m)*q(i, n)*Vq(i, n) + (q(i, m)*Vq(i, l) + q(i, l)*Vq(i, m))*(-1. + 2.*square(q(i, n))));
	}
	return ret / 6.;
}
double m1_n2_f_m1_l1_j1(int m, int n, int l, int j)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += 2.*q(i, j)*q(i, l)*q(i, n)*Vq(i, n)*square(q(i, m)) + (V(i)*q(i, j)*q(i, l) + q(i, m)*(q(i, l)*q(i, m)*Vq(i, j) + q(i, j)*q(i, m)*Vq(i, l) + q(i, j)*q(i, l)*Vq(i, m)))*(-1. + 2.*square(q(i, n)));
	}
	return ret / (2.*sqrt(2.));
}
double m2_n1_f_m2_l1(int m, int n, int l)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += q(i, n)*(Vq(i, l) + 2.*(2.*q(i, m)*(q(i, m)*Vq(i, l) + q(i, l)*Vq(i, m))*(-1. + square(q(i, m))) + V(i)*q(i, l)*(-1. + 6.*square(q(i, m))))) + q(i, l)*Vq(i, n)*square(1. - 2.*square(q(i, m)));
	}
	return ret / 6.;
}

double m2_n1_f_m1_l1_j1(int m, int n, int l, int j)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += 4.0*V(i)*q(i, j)*q(i, l)*q(i, m)*q(i, n) + (q(i, l)*q(i, m)*q(i, n)*Vq(i, j) + q(i, j)*(q(i, m)*q(i, n)*Vq(i, l) + q(i, l)*q(i, n)*Vq(i, m) + q(i, l)*q(i, m)*Vq(i, n)))*(-1.0 + 2.0*square(q(i, m)));
	}
	return ret / (2.*sqrt(2.0));
}
double m1_n1_p1_f_m1_l1_j1(int m, int n, int p, int l, int j)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += 4.*(V(i)*q(i, j)*q(i, l)*q(i, n)*q(i, p) + q(i, m)*(q(i, l)*q(i, m)*q(i, n)*q(i, p)*Vq(i, j) + q(i, j)*(q(i, m)*q(i, n)*q(i, p)*Vq(i, l) + q(i, l)*q(i, n)*q(i, p)*Vq(i, m) + q(i, l)*q(i, m)*q(i, p)*Vq(i, n) + q(i, l)*q(i, m)*q(i, n)*Vq(i, p))));
	}
	return ret / 5.;
}
double m1_n1_f_m1_n1(int m, int n)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += q(i, m)*(V(i)*q(i, m) + q(i, n)*(q(i, n)*Vq(i, m) + q(i, m)*Vq(i, n))) + V(i)*square(q(i, n));
	}
	return ret;
}
double m1_n1_f_m1_n2(int m, int n)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += q(i, n)*(V(i) + q(i, m)*Vq(i, m))*(-1. + 2.*square(q(i, n))) + square(q(i, m))*(4.*V(i)*q(i, n) + Vq(i, n)*(-1. + 2.*square(q(i, n))));
	}
	return ret / 2.;
}
double m1_n1_f_m2_n1(int m, int n)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += (V(i)*q(i, m) + q(i, n)*(q(i, n)*Vq(i, m) + q(i, m)*Vq(i, n)))*(-1. + 2.*square(q(i, m))) + 4.*V(i)*q(i, m)*square(q(i, n));
	}
	return ret / 2.;
}
double m1_n1_f_m1_n1_j1(int m, int n, int j)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += 2.*sqrt(2.)*(q(i, m)*(V(i)*q(i, j)*q(i, m) + q(i, n)*(q(i, m)*q(i, n)*Vq(i, j) + q(i, j)*q(i, n)*Vq(i, m) + q(i, j)*q(i, m)*Vq(i, n))) + V(i)*q(i, j)*square(q(i, n)));
	}
	return ret / 3.;
}
double m1_n2_f_m1_n2(int m, int n)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += square(q(i, m))*(2.*q(i, n)*Vq(i, n)*(-1. + square(q(i, n))) + V(i)*(-1. + 6.*square(q(i, n)))) + ((V(i) + q(i, m)*Vq(i, m))*square(1. - 2.*square(q(i, n))))/2.;
	}
	return ret / 2.;
}
double m1_n2_f_m2_n1(int m, int n)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += 4.*V(i)*q(i, m)*q(i, n)*(-1. + 2.*square(q(i, n))) + (-1. + 2.*square(q(i, m)))*(4.*V(i)*q(i, m)*q(i, n) + (q(i, n)*Vq(i, m) + q(i, m)*Vq(i, n))*(-1. + 2.*square(q(i, n))));
	}
	return ret / 4.;
}
double m1_n2_f_m1_n1_j1(int m, int n, int j)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += sqrt(2.)*(V(i)*q(i, j)*q(i, n)*(-1. + 2.*square(q(i, n))) + q(i, m)*(4.*V(i)*q(i, j)*q(i, m)*q(i, n) + (q(i, m)*q(i, n)*Vq(i, j) + q(i, j)*q(i, n)*Vq(i, m) + q(i, j)*q(i, m)*Vq(i, n))*(-1. + 2.*square(q(i, n)))));
	}
	return ret / 3.;
}
double m2_n1_f_m2_n1(int m, int n)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += ((V(i) + q(i, n)*Vq(i, n))*square(1. - 2.*square(q(i, m))))/2. + (2.*q(i, m)*Vq(i, m)*(-1. + square(q(i, m))) + V(i)*(-1. + 6.*square(q(i, m))))*square(q(i, n));
	}
	return ret / 2.;
}
double m2_n1_f_m1_n1_j1(int m, int n, int j)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += sqrt(2.)*((V(i)*q(i, j)*q(i, m) + q(i, n)*(q(i, m)*q(i, n)*Vq(i, j) + q(i, j)*q(i, n)*Vq(i, m) + q(i, j)*q(i, m)*Vq(i, n)))*(-1. + 2.*square(q(i, m))) + 4.*V(i)*q(i, j)*q(i, m)*square(q(i, n)));
	}
	return ret / 3.;
}
double m1_n1_p1_f_m1_n1_j1(int m, int n, int p, int j)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += q(i, m)*(V(i)*q(i, j)*q(i, m)*q(i, p) + q(i, n)*(q(i, m)*q(i, n)*q(i, p)*Vq(i, j) + q(i, j)*q(i, n)*q(i, p)*Vq(i, m) + q(i, j)*q(i, m)*q(i, p)*Vq(i, n) + q(i, j)*q(i, m)*q(i, n)*Vq(i, p))) + V(i)*q(i, j)*q(i, p)*square(q(i, n));
	}
	return ret;
}
double m1_n1_p1_f_m1_n1_p1(int m, int n, int p)
{
    double ret = 0;
	for (int i=0; i < V.n_rows; i++) {
		ret += 4.*((V(i) + q(i, p)*Vq(i, p))*square(q(i, m))*square(q(i, n)) + (V(i) + q(i, n)*Vq(i, n))*square(q(i, m))*square(q(i, p)) + (V(i) + q(i, m)*Vq(i, m))*square(q(i, n))*square(q(i, p)));
	}
	return ret / 3.;
}
