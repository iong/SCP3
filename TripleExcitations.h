double f_k1(int k)
{
    double *Vqk = (double *) __builtin_assume_aligned(Vq.colptr(k), 32);
    double * qk = (double *) __builtin_assume_aligned( q.colptr(k), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        ret += Vqk[i];
    }
    return ret / sqrt(2.0);
}
double f_k2(int k)
{
    double *Vqk = (double *) __builtin_assume_aligned(Vq.colptr(k), 32);
    double * qk = (double *) __builtin_assume_aligned( q.colptr(k), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        ret += qk[i]*Vqk[i];
    }
    return ret / sqrt(2.);
}
double f_k3(int k)
{
    double *Vqk = (double *) __builtin_assume_aligned(Vq.colptr(k), 32);
    double * qk = (double *) __builtin_assume_aligned( q.colptr(k), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        double _sqr_qki = qk[i]*qk[i];
        ret += 4.*V[i]*qk[i] + Vqk[i]*(-3. + 2.*_sqr_qki);
    }
    return ret / (2.*sqrt(3.));
}
double f_k1_l1(int k, int l)
{
    double *Vqk = (double *) __builtin_assume_aligned(Vq.colptr(k), 32);
    double * qk = (double *) __builtin_assume_aligned( q.colptr(k), 32);
    double *Vql = (double *) __builtin_assume_aligned(Vq.colptr(l), 32);
    double * ql = (double *) __builtin_assume_aligned( q.colptr(l), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        ret += ql[i]*Vqk[i] + qk[i]*Vql[i];
    }
    return ret / 2.;
}
double f_k1_l2(int k, int l)
{
    double *Vqk = (double *) __builtin_assume_aligned(Vq.colptr(k), 32);
    double * qk = (double *) __builtin_assume_aligned( q.colptr(k), 32);
    double *Vql = (double *) __builtin_assume_aligned(Vq.colptr(l), 32);
    double * ql = (double *) __builtin_assume_aligned( q.colptr(l), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        double _sqr_qli = ql[i]*ql[i];
        ret += qk[i]*ql[i]*Vql[i] + Vqk[i]*(-0.5 + _sqr_qli);
    }
    return ret / 2.;
}
double f_k2_l1(int k, int l)
{
    double *Vqk = (double *) __builtin_assume_aligned(Vq.colptr(k), 32);
    double * qk = (double *) __builtin_assume_aligned( q.colptr(k), 32);
    double *Vql = (double *) __builtin_assume_aligned(Vq.colptr(l), 32);
    double * ql = (double *) __builtin_assume_aligned( q.colptr(l), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        double _sqr_qki = qk[i]*qk[i];
        ret += qk[i]*ql[i]*Vqk[i] + Vql[i]*(-0.5 + _sqr_qki);
    }
    return ret / 2.;
}
double f_k1_l1_j1(int k, int l, int j)
{
    double *Vqk = (double *) __builtin_assume_aligned(Vq.colptr(k), 32);
    double * qk = (double *) __builtin_assume_aligned( q.colptr(k), 32);
    double *Vql = (double *) __builtin_assume_aligned(Vq.colptr(l), 32);
    double * ql = (double *) __builtin_assume_aligned( q.colptr(l), 32);
    double *Vqj = (double *) __builtin_assume_aligned(Vq.colptr(j), 32);
    double * qj = (double *) __builtin_assume_aligned( q.colptr(j), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        ret += sqrt(2.)*(qk[i]*ql[i]*Vqj[i] + qj[i]*ql[i]*Vqk[i] + qj[i]*qk[i]*Vql[i]);
    }
    return ret / 3.;
}

double m1_f_k1(int m, int k)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double *Vqk = (double *) __builtin_assume_aligned(Vq.colptr(k), 32);
    double * qk = (double *) __builtin_assume_aligned( q.colptr(k), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        ret += qm[i]*Vqk[i] + qk[i]*Vqm[i];
    }
    return ret / 2.;
}

double m1_f_k2(int m, int k)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double *Vqk = (double *) __builtin_assume_aligned(Vq.colptr(k), 32);
    double * qk = (double *) __builtin_assume_aligned( q.colptr(k), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        double _sqr_qki = qk[i]*qk[i];
        ret += qk[i]*qm[i]*Vqk[i] + Vqm[i]*(-0.5 + _sqr_qki);
    }
    return ret / 2.;
}


double m1_f_k3(int m, int k)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double *Vqk = (double *) __builtin_assume_aligned(Vq.colptr(k), 32);
    double * qk = (double *) __builtin_assume_aligned( q.colptr(k), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        double _sqr_qki = qk[i]*qk[i];
        ret += 4.*V[i]*qk[i]*qm[i] + (qm[i]*Vqk[i] + qk[i]*Vqm[i])*(-3. + 2.*_sqr_qki);
    }
    return ret / (2.*sqrt(6.));
}


double m1_f_k1_l1(int m, int k, int l)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double *Vqk = (double *) __builtin_assume_aligned(Vq.colptr(k), 32);
    double * qk = (double *) __builtin_assume_aligned( q.colptr(k), 32);
    double *Vql = (double *) __builtin_assume_aligned(Vq.colptr(l), 32);
    double * ql = (double *) __builtin_assume_aligned( q.colptr(l), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        ret += sqrt(2.)*(ql[i]*qm[i]*Vqk[i] + qk[i]*qm[i]*Vql[i] + qk[i]*ql[i]*Vqm[i]);
    }
    return ret / 3.;
}


double m1_f_k1_l2(int m, int k, int l)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double *Vqk = (double *) __builtin_assume_aligned(Vq.colptr(k), 32);
    double * qk = (double *) __builtin_assume_aligned( q.colptr(k), 32);
    double *Vql = (double *) __builtin_assume_aligned(Vq.colptr(l), 32);
    double * ql = (double *) __builtin_assume_aligned( q.colptr(l), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        double _sqr_qli = ql[i]*ql[i];
        ret += 2.*qk[i]*ql[i]*qm[i]*Vql[i] + qm[i]*Vqk[i]*(-1. + 2.*_sqr_qli) + qk[i]*Vqm[i]*(-1. + 2.*_sqr_qli);
    }
    return ret / (3.*sqrt(2.));
}
double m1_f_k2_l1(int m, int k, int l)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double *Vqk = (double *) __builtin_assume_aligned(Vq.colptr(k), 32);
    double * qk = (double *) __builtin_assume_aligned( q.colptr(k), 32);
    double *Vql = (double *) __builtin_assume_aligned(Vq.colptr(l), 32);
    double * ql = (double *) __builtin_assume_aligned( q.colptr(l), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        double _sqr_qki = qk[i]*qk[i];
        ret += 2.*qk[i]*ql[i]*qm[i]*Vqk[i] + (qm[i]*Vql[i] + ql[i]*Vqm[i])*(-1. + 2.*_sqr_qki);
    }
    return ret / (3.*sqrt(2.));
}
double m1_f_k1_l1_j1(int m, int k, int l, int j)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double *Vqk = (double *) __builtin_assume_aligned(Vq.colptr(k), 32);
    double * qk = (double *) __builtin_assume_aligned( q.colptr(k), 32);
    double *Vql = (double *) __builtin_assume_aligned(Vq.colptr(l), 32);
    double * ql = (double *) __builtin_assume_aligned( q.colptr(l), 32);
    double *Vqj = (double *) __builtin_assume_aligned(Vq.colptr(j), 32);
    double * qj = (double *) __builtin_assume_aligned( q.colptr(j), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        ret += qk[i]*ql[i]*qm[i]*Vqj[i] + qj[i]*(ql[i]*qm[i]*Vqk[i] + qk[i]*qm[i]*Vql[i] + qk[i]*ql[i]*Vqm[i]);
    }
    return ret / 2.;
}
double m2_f_k2(int m, int k)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double *Vqk = (double *) __builtin_assume_aligned(Vq.colptr(k), 32);
    double * qk = (double *) __builtin_assume_aligned( q.colptr(k), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        double _sqr_qki = qk[i]*qk[i];
        double _sqr_qmi = qm[i]*qm[i];
        ret += qm[i]*Vqm[i]*(-1. + 2.*_sqr_qki) + qk[i]*Vqk[i]*(-1. + 2.*_sqr_qmi);
    }
    return ret / 4.;
}
double m2_f_k3(int m, int k)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double *Vqk = (double *) __builtin_assume_aligned(Vq.colptr(k), 32);
    double * qk = (double *) __builtin_assume_aligned( q.colptr(k), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        double _sqr_qki = qk[i]*qk[i];
        double _sqr_qmi = qm[i]*qm[i];
        ret += 4.*V[i]*qk[i]*(-1. + 2.*_sqr_qmi) + (-3. + 2.*_sqr_qki)*(2.*qk[i]*qm[i]*Vqm[i] + Vqk[i]*(-1. + 2.*_sqr_qmi));
    }
    return ret / (4.*sqrt(6.));
}
double m2_f_k1_l1(int m, int k, int l)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double *Vqk = (double *) __builtin_assume_aligned(Vq.colptr(k), 32);
    double * qk = (double *) __builtin_assume_aligned( q.colptr(k), 32);
    double *Vql = (double *) __builtin_assume_aligned(Vq.colptr(l), 32);
    double * ql = (double *) __builtin_assume_aligned( q.colptr(l), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        double _sqr_qmi = qm[i]*qm[i];
        ret += 2.*qk[i]*ql[i]*qm[i]*Vqm[i] + (ql[i]*Vqk[i] + qk[i]*Vql[i])*(-1. + 2.*_sqr_qmi);
    }
    return ret / (3.*sqrt(2.));
}
double m2_f_k1_l2(int m, int k, int l)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double *Vqk = (double *) __builtin_assume_aligned(Vq.colptr(k), 32);
    double * qk = (double *) __builtin_assume_aligned( q.colptr(k), 32);
    double *Vql = (double *) __builtin_assume_aligned(Vq.colptr(l), 32);
    double * ql = (double *) __builtin_assume_aligned( q.colptr(l), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        double _sqr_qmi = qm[i]*qm[i];
        double _sqr_qli = ql[i]*ql[i];
        ret += 2.*qk[i]*qm[i]*Vqm[i]*(-1. + 2.*_sqr_qli) + (2.*qk[i]*ql[i]*Vql[i] + Vqk[i]*(-1. + 2.*_sqr_qli))*(-1. + 2.*_sqr_qmi);
    }
    return ret / (6.*sqrt(2.));
}
double m2_f_k2_l1(int m, int k, int l)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double *Vqk = (double *) __builtin_assume_aligned(Vq.colptr(k), 32);
    double * qk = (double *) __builtin_assume_aligned( q.colptr(k), 32);
    double *Vql = (double *) __builtin_assume_aligned(Vq.colptr(l), 32);
    double * ql = (double *) __builtin_assume_aligned( q.colptr(l), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        double _sqr_qki = qk[i]*qk[i];
        double _sqr_qmi = qm[i]*qm[i];
        ret += 2.*qk[i]*ql[i]*Vqk[i]*(-1. + 2.*_sqr_qmi) + (-1. + 2.*_sqr_qki)*(2.*ql[i]*qm[i]*Vqm[i] + Vql[i]*(-1. + 2.*_sqr_qmi));
    }
    return ret / (6.*sqrt(2.));
}
double m2_f_k1_l1_j1(int m, int k, int l, int j)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double *Vqk = (double *) __builtin_assume_aligned(Vq.colptr(k), 32);
    double * qk = (double *) __builtin_assume_aligned( q.colptr(k), 32);
    double *Vql = (double *) __builtin_assume_aligned(Vq.colptr(l), 32);
    double * ql = (double *) __builtin_assume_aligned( q.colptr(l), 32);
    double *Vqj = (double *) __builtin_assume_aligned(Vq.colptr(j), 32);
    double * qj = (double *) __builtin_assume_aligned( q.colptr(j), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        double _sqr_qmi = qm[i]*qm[i];
        ret += 2.*qj[i]*qk[i]*ql[i]*qm[i]*Vqm[i] + (qk[i]*ql[i]*Vqj[i] + qj[i]*ql[i]*Vqk[i] + qj[i]*qk[i]*Vql[i])*(-1. + 2.*_sqr_qmi);
    }
    return ret / 4.;
}
double m3_f_k3(int m, int k)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double *Vqk = (double *) __builtin_assume_aligned(Vq.colptr(k), 32);
    double * qk = (double *) __builtin_assume_aligned( q.colptr(k), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        double _sqr_qki = qk[i]*qk[i];
        double _sqr_qmi = qm[i]*qm[i];
        ret += 4.*V[i]*qk[i]*qm[i]*(-3. + 2.*_sqr_qmi) + (-3. + 2.*_sqr_qki)*(4.*V[i]*qk[i]*qm[i] + (qm[i]*Vqk[i] + qk[i]*Vqm[i])*(-3. + 2.*_sqr_qmi));
    }
    return ret / 12.;
}
double m3_f_k1_l1(int m, int k, int l)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double *Vqk = (double *) __builtin_assume_aligned(Vq.colptr(k), 32);
    double * qk = (double *) __builtin_assume_aligned( q.colptr(k), 32);
    double *Vql = (double *) __builtin_assume_aligned(Vq.colptr(l), 32);
    double * ql = (double *) __builtin_assume_aligned( q.colptr(l), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        double _sqr_qmi = qm[i]*qm[i];
        ret += 4.*V[i]*qk[i]*ql[i]*qm[i] + (ql[i]*qm[i]*Vqk[i] + qk[i]*qm[i]*Vql[i] + qk[i]*ql[i]*Vqm[i])*(-3. + 2.*_sqr_qmi);
    }
    return ret / (3.*sqrt(3.));
}
double m3_f_k1_l2(int m, int k, int l)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double *Vqk = (double *) __builtin_assume_aligned(Vq.colptr(k), 32);
    double * qk = (double *) __builtin_assume_aligned( q.colptr(k), 32);
    double *Vql = (double *) __builtin_assume_aligned(Vq.colptr(l), 32);
    double * ql = (double *) __builtin_assume_aligned( q.colptr(l), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        double _sqr_qmi = qm[i]*qm[i];
        double _sqr_qli = ql[i]*ql[i];
        ret += 4.*V[i]*qk[i]*qm[i]*(-1. + 2.*_sqr_qli) + (2.*qk[i]*ql[i]*qm[i]*Vql[i] + qm[i]*Vqk[i]*(-1. + 2.*_sqr_qli) + qk[i]*Vqm[i]*(-1. + 2.*_sqr_qli))*(-3. + 2.*_sqr_qmi);
    }
    return ret / (6.*sqrt(3.));
}
double m3_f_k2_l1(int m, int k, int l)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double *Vqk = (double *) __builtin_assume_aligned(Vq.colptr(k), 32);
    double * qk = (double *) __builtin_assume_aligned( q.colptr(k), 32);
    double *Vql = (double *) __builtin_assume_aligned(Vq.colptr(l), 32);
    double * ql = (double *) __builtin_assume_aligned( q.colptr(l), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        double _sqr_qki = qk[i]*qk[i];
        double _sqr_qmi = qm[i]*qm[i];
        ret += 4.*V[i]*ql[i]*qm[i]*(-1. + 2.*_sqr_qki) + (2.*qk[i]*ql[i]*qm[i]*Vqk[i] + (qm[i]*Vql[i] + ql[i]*Vqm[i])*(-1. + 2.*_sqr_qki))*(-3. + 2.*_sqr_qmi);
    }
    return ret / (6.*sqrt(3.));
}
double m3_f_k1_l1_j1(int m, int k, int l, int j)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double *Vqk = (double *) __builtin_assume_aligned(Vq.colptr(k), 32);
    double * qk = (double *) __builtin_assume_aligned( q.colptr(k), 32);
    double *Vql = (double *) __builtin_assume_aligned(Vq.colptr(l), 32);
    double * ql = (double *) __builtin_assume_aligned( q.colptr(l), 32);
    double *Vqj = (double *) __builtin_assume_aligned(Vq.colptr(j), 32);
    double * qj = (double *) __builtin_assume_aligned( q.colptr(j), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        double _sqr_qmi = qm[i]*qm[i];
        ret += 4.*V[i]*qj[i]*qk[i]*ql[i]*qm[i] + (qk[i]*ql[i]*qm[i]*Vqj[i] + qj[i]*(ql[i]*qm[i]*Vqk[i] + qk[i]*qm[i]*Vql[i] + qk[i]*ql[i]*Vqm[i]))*(-3. + 2.*_sqr_qmi);
    }
    return ret / (2.*sqrt(6.));
}
double m1_n1_f_k1_l1(int m, int n, int k, int l)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double *Vqn = (double *) __builtin_assume_aligned(Vq.colptr(n), 32);
    double * qn = (double *) __builtin_assume_aligned( q.colptr(n), 32);
    double *Vqk = (double *) __builtin_assume_aligned(Vq.colptr(k), 32);
    double * qk = (double *) __builtin_assume_aligned( q.colptr(k), 32);
    double *Vql = (double *) __builtin_assume_aligned(Vq.colptr(l), 32);
    double * ql = (double *) __builtin_assume_aligned( q.colptr(l), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        ret += ql[i]*qm[i]*qn[i]*Vqk[i] + qk[i]*(qm[i]*qn[i]*Vql[i] + ql[i]*qn[i]*Vqm[i] + ql[i]*qm[i]*Vqn[i]);
    }
    return ret / 2.;
}
double m1_n1_f_k1_l2(int m, int n, int k, int l)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double *Vqn = (double *) __builtin_assume_aligned(Vq.colptr(n), 32);
    double * qn = (double *) __builtin_assume_aligned( q.colptr(n), 32);
    double *Vqk = (double *) __builtin_assume_aligned(Vq.colptr(k), 32);
    double * qk = (double *) __builtin_assume_aligned( q.colptr(k), 32);
    double *Vql = (double *) __builtin_assume_aligned(Vq.colptr(l), 32);
    double * ql = (double *) __builtin_assume_aligned( q.colptr(l), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        double _sqr_qli = ql[i]*ql[i];
        ret += 2.*qk[i]*ql[i]*qm[i]*qn[i]*Vql[i] + qm[i]*qn[i]*Vqk[i]*(-1. + 2.*_sqr_qli) + qk[i]*(qn[i]*Vqm[i] + qm[i]*Vqn[i])*(-1. + 2.*_sqr_qli);
    }
    return ret / 4.;
}
double m1_n1_f_k2_l1(int m, int n, int k, int l)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double *Vqn = (double *) __builtin_assume_aligned(Vq.colptr(n), 32);
    double * qn = (double *) __builtin_assume_aligned( q.colptr(n), 32);
    double *Vqk = (double *) __builtin_assume_aligned(Vq.colptr(k), 32);
    double * qk = (double *) __builtin_assume_aligned( q.colptr(k), 32);
    double *Vql = (double *) __builtin_assume_aligned(Vq.colptr(l), 32);
    double * ql = (double *) __builtin_assume_aligned( q.colptr(l), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        double _sqr_qki = qk[i]*qk[i];
        ret += 2.*qk[i]*ql[i]*qm[i]*qn[i]*Vqk[i] + (qm[i]*qn[i]*Vql[i] + ql[i]*qn[i]*Vqm[i] + ql[i]*qm[i]*Vqn[i])*(-1. + 2.*_sqr_qki);
    }
    return ret / 4.;
}
double m1_n1_f_k1_l1_j1(int m, int n, int k, int l, int j)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double *Vqn = (double *) __builtin_assume_aligned(Vq.colptr(n), 32);
    double * qn = (double *) __builtin_assume_aligned( q.colptr(n), 32);
    double *Vqk = (double *) __builtin_assume_aligned(Vq.colptr(k), 32);
    double * qk = (double *) __builtin_assume_aligned( q.colptr(k), 32);
    double *Vql = (double *) __builtin_assume_aligned(Vq.colptr(l), 32);
    double * ql = (double *) __builtin_assume_aligned( q.colptr(l), 32);
    double *Vqj = (double *) __builtin_assume_aligned(Vq.colptr(j), 32);
    double * qj = (double *) __builtin_assume_aligned( q.colptr(j), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        ret += 2.*sqrt(2.)*(qk[i]*ql[i]*qm[i]*qn[i]*Vqj[i] + qj[i]*(ql[i]*qm[i]*qn[i]*Vqk[i] + qk[i]*(qm[i]*qn[i]*Vql[i] + ql[i]*qn[i]*Vqm[i] + ql[i]*qm[i]*Vqn[i])));
    }
    return ret / 5.;
}
double m1_n2_f_k1_l2(int m, int n, int k, int l)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double *Vqn = (double *) __builtin_assume_aligned(Vq.colptr(n), 32);
    double * qn = (double *) __builtin_assume_aligned( q.colptr(n), 32);
    double *Vqk = (double *) __builtin_assume_aligned(Vq.colptr(k), 32);
    double * qk = (double *) __builtin_assume_aligned( q.colptr(k), 32);
    double *Vql = (double *) __builtin_assume_aligned(Vq.colptr(l), 32);
    double * ql = (double *) __builtin_assume_aligned( q.colptr(l), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        double _sqr_qli = ql[i]*ql[i];
        double _sqr_qni = qn[i]*qn[i];
        ret += 2.*qk[i]*qm[i]*qn[i]*Vqn[i]*(-1. + 2.*_sqr_qli) + (2.*qk[i]*ql[i]*qm[i]*Vql[i] + qm[i]*Vqk[i]*(-1. + 2.*_sqr_qli) + qk[i]*Vqm[i]*(-1. + 2.*_sqr_qli))*(-1. + 2.*_sqr_qni);
    }
    return ret / 8.;
}
double m1_n2_f_k2_l1(int m, int n, int k, int l)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double *Vqn = (double *) __builtin_assume_aligned(Vq.colptr(n), 32);
    double * qn = (double *) __builtin_assume_aligned( q.colptr(n), 32);
    double *Vqk = (double *) __builtin_assume_aligned(Vq.colptr(k), 32);
    double * qk = (double *) __builtin_assume_aligned( q.colptr(k), 32);
    double *Vql = (double *) __builtin_assume_aligned(Vq.colptr(l), 32);
    double * ql = (double *) __builtin_assume_aligned( q.colptr(l), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        double _sqr_qki = qk[i]*qk[i];
        double _sqr_qni = qn[i]*qn[i];
        ret += 2.*qk[i]*ql[i]*qm[i]*Vqk[i]*(-1. + 2.*_sqr_qni) + (-1. + 2.*_sqr_qki)*(2.*ql[i]*qm[i]*qn[i]*Vqn[i] + (qm[i]*Vql[i] + ql[i]*Vqm[i])*(-1. + 2.*_sqr_qni));
    }
    return ret / 8.;
}
double m1_n2_f_k1_l1_j1(int m, int n, int k, int l, int j)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double *Vqn = (double *) __builtin_assume_aligned(Vq.colptr(n), 32);
    double * qn = (double *) __builtin_assume_aligned( q.colptr(n), 32);
    double *Vqk = (double *) __builtin_assume_aligned(Vq.colptr(k), 32);
    double * qk = (double *) __builtin_assume_aligned( q.colptr(k), 32);
    double *Vql = (double *) __builtin_assume_aligned(Vq.colptr(l), 32);
    double * ql = (double *) __builtin_assume_aligned( q.colptr(l), 32);
    double *Vqj = (double *) __builtin_assume_aligned(Vq.colptr(j), 32);
    double * qj = (double *) __builtin_assume_aligned( q.colptr(j), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        double _sqr_qni = qn[i]*qn[i];
        ret += sqrt(2.)*(2.*qj[i]*qk[i]*ql[i]*qm[i]*qn[i]*Vqn[i] + (qk[i]*ql[i]*qm[i]*Vqj[i] + qj[i]*(ql[i]*qm[i]*Vqk[i] + qk[i]*qm[i]*Vql[i] + qk[i]*ql[i]*Vqm[i]))*(-1. + 2.*_sqr_qni));
    }
    return ret / 5.;
}
double m2_n1_f_k2_l1(int m, int n, int k, int l)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double *Vqn = (double *) __builtin_assume_aligned(Vq.colptr(n), 32);
    double * qn = (double *) __builtin_assume_aligned( q.colptr(n), 32);
    double *Vqk = (double *) __builtin_assume_aligned(Vq.colptr(k), 32);
    double * qk = (double *) __builtin_assume_aligned( q.colptr(k), 32);
    double *Vql = (double *) __builtin_assume_aligned(Vq.colptr(l), 32);
    double * ql = (double *) __builtin_assume_aligned( q.colptr(l), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        double _sqr_qki = qk[i]*qk[i];
        double _sqr_qmi = qm[i]*qm[i];
        ret += 2.*qk[i]*ql[i]*qn[i]*Vqk[i]*(-1. + 2.*_sqr_qmi) + (-1. + 2.*_sqr_qki)*(2.*ql[i]*qm[i]*qn[i]*Vqm[i] + qn[i]*Vql[i]*(-1. + 2.*_sqr_qmi) + ql[i]*Vqn[i]*(-1. + 2.*_sqr_qmi));
    }
    return ret / 8.;
}
double m2_n1_f_k1_l1_j1(int m, int n, int k, int l, int j)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double *Vqn = (double *) __builtin_assume_aligned(Vq.colptr(n), 32);
    double * qn = (double *) __builtin_assume_aligned( q.colptr(n), 32);
    double *Vqk = (double *) __builtin_assume_aligned(Vq.colptr(k), 32);
    double * qk = (double *) __builtin_assume_aligned( q.colptr(k), 32);
    double *Vql = (double *) __builtin_assume_aligned(Vq.colptr(l), 32);
    double * ql = (double *) __builtin_assume_aligned( q.colptr(l), 32);
    double *Vqj = (double *) __builtin_assume_aligned(Vq.colptr(j), 32);
    double * qj = (double *) __builtin_assume_aligned( q.colptr(j), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        double _sqr_qmi = qm[i]*qm[i];
        ret += sqrt(2.)*(2.*qj[i]*qk[i]*ql[i]*qm[i]*qn[i]*Vqm[i] + qn[i]*(qk[i]*ql[i]*Vqj[i] + qj[i]*ql[i]*Vqk[i] + qj[i]*qk[i]*Vql[i])*(-1. + 2.*_sqr_qmi) + qj[i]*qk[i]*ql[i]*Vqn[i]*(-1. + 2.*_sqr_qmi));
    }
    return ret / 5.;
}

#ifdef __FMA4__
double m1_n1_p1_f_k1_l1_j1(int m, int n, int p, int k, int l, int j)
{
    __m256d vret  = _mm256_setzero_pd();

    double *Vqm = Vq.colptr(m);
    double *Vqn = Vq.colptr(n);
    double *Vqp = Vq.colptr(p);
    double *Vqk = Vq.colptr(k);
    double *Vql = Vq.colptr(l);
    double *Vqj = Vq.colptr(j);

    double *qm = Vq.colptr(m);
    double *qn = Vq.colptr(n);
    double *qp = Vq.colptr(p);
    double *qk = Vq.colptr(k);
    double *ql = Vq.colptr(l);
    double *qj = Vq.colptr(j);

    for (int i=0; i < V.n_rows; i += 4) {
        __m256d Vs = _mm256_load_pd(Vqm + i);
        __m256d qs = _mm256_load_pd( qm + i);

        __m256d Vq4 = _mm256_load_pd(Vqn + i);
        __m256d  q4 = _mm256_load_pd( qn + i);

        Vs = _mm256_macc_pd(Vs, q4, _mm256_mul_pd(qs, Vq4));
        qs = _mm256_mul_pd(qs, q4);

        Vq4 = _mm256_load_pd(Vqp + i);
         q4 = _mm256_load_pd( qp + i);
        Vs = _mm256_macc_pd(Vs, q4, _mm256_mul_pd(qs, Vq4));
        qs = _mm256_mul_pd(qs, q4);

        Vq4 = _mm256_load_pd(Vqk + i);
         q4 = _mm256_load_pd( qk + i);
        Vs = _mm256_macc_pd(Vs, q4, _mm256_mul_pd(qs, Vq4));
        qs = _mm256_mul_pd(qs, q4);

        Vq4 = _mm256_load_pd(Vql + i);
         q4 = _mm256_load_pd( ql + i);
        Vs = _mm256_macc_pd(Vs, q4, _mm256_mul_pd(qs, Vq4));
        qs = _mm256_mul_pd(qs, q4);

        Vq4 = _mm256_load_pd(Vqj + i);
         q4 = _mm256_load_pd( qj + i);
        Vs = _mm256_macc_pd(Vs, q4, _mm256_mul_pd(qs, Vq4));

        vret = _mm256_add_pd(vret, Vs);
    }

    double ret[4] __attribute__((aligned(32)));
    _mm256_store_pd(ret, vret);
        
    return 2.0 / 3. * (ret[0] + ret[1] + ret[2] + ret[3]);
}

#else

double m1_n1_p1_f_k1_l1_j1(int m, int n, int p, int k, int l, int j)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double *Vqn = (double *) __builtin_assume_aligned(Vq.colptr(n), 32);
    double * qn = (double *) __builtin_assume_aligned( q.colptr(n), 32);
    double *Vqp = (double *) __builtin_assume_aligned(Vq.colptr(p), 32);
    double * qp = (double *) __builtin_assume_aligned( q.colptr(p), 32);
    double *Vqk = (double *) __builtin_assume_aligned(Vq.colptr(k), 32);
    double * qk = (double *) __builtin_assume_aligned( q.colptr(k), 32);
    double *Vql = (double *) __builtin_assume_aligned(Vq.colptr(l), 32);
    double * ql = (double *) __builtin_assume_aligned( q.colptr(l), 32);
    double *Vqj = (double *) __builtin_assume_aligned(Vq.colptr(j), 32);
    double * qj = (double *) __builtin_assume_aligned( q.colptr(j), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        ret += 2.*(qk[i]*ql[i]*qm[i]*qn[i]*qp[i]*Vqj[i] + qj[i]*ql[i]*qm[i]*qn[i]*qp[i]*Vqk[i] + qj[i]*qk[i]*(qm[i]*qn[i]*qp[i]*Vql[i] + ql[i]*qn[i]*qp[i]*Vqm[i] + ql[i]*qm[i]*qp[i]*Vqn[i] + ql[i]*qm[i]*qn[i]*Vqp[i]));
    }
    return ret / 3.;
}
#endif
double m1_f_m1(int m)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        ret += V[i] + qm[i]*Vqm[i];
    }
    return ret;
}
double m1_f_m2(int m)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        double _sqr_qmi = qm[i]*qm[i];
        ret += 2.*V[i]*qm[i] + Vqm[i]*(-0.5 + _sqr_qmi);
    }
    return ret;
}
double m1_f_m3(int m)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        double _sqr_qmi = qm[i]*qm[i];
        ret += qm[i]*Vqm[i]*(-3. + 2.*_sqr_qmi) + V[i]*(-3. + 6.*_sqr_qmi);
    }
    return ret / sqrt(6.);
}
double m1_f_m1_l1(int m, int l)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double *Vql = (double *) __builtin_assume_aligned(Vq.colptr(l), 32);
    double * ql = (double *) __builtin_assume_aligned( q.colptr(l), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        ret += V[i]*ql[i] + qm[i]*(qm[i]*Vql[i] + ql[i]*Vqm[i]);
    }
    return ret / sqrt(2.);
}
double m1_f_m1_l2(int m, int l)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double *Vql = (double *) __builtin_assume_aligned(Vq.colptr(l), 32);
    double * ql = (double *) __builtin_assume_aligned( q.colptr(l), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        double _sqr_qmi = qm[i]*qm[i];
        double _sqr_qli = ql[i]*ql[i];
        ret += (V[i] + qm[i]*Vqm[i])*(-1. + 2.*_sqr_qli) + 2.*ql[i]*Vql[i]*_sqr_qmi;
    }
    return ret / (2.*sqrt(2.));
}
double m1_f_m2_l1(int m, int l)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double *Vql = (double *) __builtin_assume_aligned(Vq.colptr(l), 32);
    double * ql = (double *) __builtin_assume_aligned( q.colptr(l), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        double _sqr_qmi = qm[i]*qm[i];
        ret += 4.*V[i]*ql[i]*qm[i] + (qm[i]*Vql[i] + ql[i]*Vqm[i])*(-1. + 2.*_sqr_qmi);
    }
    return ret / (2.*sqrt(2.));
}
double m1_f_m1_l1_j1(int m, int l, int j)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double *Vql = (double *) __builtin_assume_aligned(Vq.colptr(l), 32);
    double * ql = (double *) __builtin_assume_aligned( q.colptr(l), 32);
    double *Vqj = (double *) __builtin_assume_aligned(Vq.colptr(j), 32);
    double * qj = (double *) __builtin_assume_aligned( q.colptr(j), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        ret += 2.*(V[i]*qj[i]*ql[i] + qm[i]*(ql[i]*qm[i]*Vqj[i] + qj[i]*qm[i]*Vql[i] + qj[i]*ql[i]*Vqm[i]));
    }
    return ret / 3.;
}
double m2_f_m2(int m)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        double _sqr_qmi = qm[i]*qm[i];
        ret += qm[i]*Vqm[i]*(-1. + _sqr_qmi) + V[i]*(-0.5 + 3.*_sqr_qmi);
    }
    return ret;
}

double m2_f_m3(int m)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        double _sqr_qmi = qm[i]*qm[i];
        ret += Vqm[i]*(3. + 4.*square(_sqr_qmi) - 8.*_sqr_qmi) + 16.*V[i]*qm[i]*(-1. + _sqr_qmi);
    }
    return ret / (2.*sqrt(6.));
}
double m2_f_m1_l1(int m, int l)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double *Vql = (double *) __builtin_assume_aligned(Vq.colptr(l), 32);
    double * ql = (double *) __builtin_assume_aligned( q.colptr(l), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        double _sqr_qmi = qm[i]*qm[i];
        ret += 4.*V[i]*ql[i]*qm[i] + (qm[i]*Vql[i] + ql[i]*Vqm[i])*(-1. + 2.*_sqr_qmi);
    }
    return ret / (2.*sqrt(2.));
}
double m2_f_m1_l2(int m, int l)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double *Vql = (double *) __builtin_assume_aligned(Vq.colptr(l), 32);
    double * ql = (double *) __builtin_assume_aligned( q.colptr(l), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        double _sqr_qmi = qm[i]*qm[i];
        double _sqr_qli = ql[i]*ql[i];
        ret += 4.*V[i]*qm[i]*(-1. + 2.*_sqr_qli) + (2.*ql[i]*qm[i]*Vql[i] + Vqm[i]*(-1. + 2.*_sqr_qli))*(-1. + 2.*_sqr_qmi);
    }
    return ret / (4.*sqrt(2.));
}
double m2_f_m2_l1(int m, int l)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double *Vql = (double *) __builtin_assume_aligned(Vq.colptr(l), 32);
    double * ql = (double *) __builtin_assume_aligned( q.colptr(l), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        double _sqr_qmi = qm[i]*qm[i];
        ret += Vql[i] + 2.*(2.*qm[i]*(qm[i]*Vql[i] + ql[i]*Vqm[i])*(-1. + _sqr_qmi) + V[i]*ql[i]*(-1. + 6.*_sqr_qmi));
    }
    return ret / (4.*sqrt(2.));
}
double m2_f_m1_l1_j1(int m, int l, int j)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double *Vql = (double *) __builtin_assume_aligned(Vq.colptr(l), 32);
    double * ql = (double *) __builtin_assume_aligned( q.colptr(l), 32);
    double *Vqj = (double *) __builtin_assume_aligned(Vq.colptr(j), 32);
    double * qj = (double *) __builtin_assume_aligned( q.colptr(j), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        double _sqr_qmi = qm[i]*qm[i];
        ret += 4.*V[i]*qj[i]*ql[i]*qm[i] + (ql[i]*qm[i]*Vqj[i] + qj[i]*qm[i]*Vql[i] + qj[i]*ql[i]*Vqm[i])*(-1. + 2.*_sqr_qmi);
    }
    return ret / 3.;
}
double m3_f_m3(int m)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        double _sqr_qmi = qm[i]*qm[i];
        ret += (-3. + 2.*_sqr_qmi)*(qm[i]*Vqm[i]*(-3. + 2.*_sqr_qmi) + V[i]*(-3. + 1.0*_sqr_qmi));
    }
    return ret / 6.;
}
double m3_f_m1_l1(int m, int l)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double *Vql = (double *) __builtin_assume_aligned(Vq.colptr(l), 32);
    double * ql = (double *) __builtin_assume_aligned( q.colptr(l), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        double _sqr_qmi = qm[i]*qm[i];
        ret += qm[i]*(qm[i]*Vql[i] + ql[i]*Vqm[i])*(-3. + 2.*_sqr_qmi) + 3.*V[i]*ql[i]*(-1. + 2.*_sqr_qmi);
    }
    return ret / (2.*sqrt(3.));
}
double m3_f_m1_l2(int m, int l)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double *Vql = (double *) __builtin_assume_aligned(Vq.colptr(l), 32);
    double * ql = (double *) __builtin_assume_aligned( q.colptr(l), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        double _sqr_qmi = qm[i]*qm[i];
        double _sqr_qli = ql[i]*ql[i];
        ret += qm[i]*(2.*ql[i]*qm[i]*Vql[i] + Vqm[i]*(-1. + 2.*_sqr_qli))*(-3. + 2.*_sqr_qmi) + 3.*V[i]*(-1. + 2.*_sqr_qli)*(-1. + 2.*_sqr_qmi);
    }
    return ret / (4.*sqrt(3.));
}
double m3_f_m2_l1(int m, int l)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double *Vql = (double *) __builtin_assume_aligned(Vq.colptr(l), 32);
    double * ql = (double *) __builtin_assume_aligned( q.colptr(l), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        double _sqr_qmi = qm[i]*qm[i];
        ret += (qm[i]*Vql[i] + ql[i]*Vqm[i])*(3. + 4.*square(_sqr_qmi) - 8.*_sqr_qmi) + 16.*V[i]*ql[i]*qm[i]*(-1. + _sqr_qmi);
    }
    return ret / (4.*sqrt(3.));
}
double m3_f_m1_l1_j1(int m, int l, int j)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double *Vql = (double *) __builtin_assume_aligned(Vq.colptr(l), 32);
    double * ql = (double *) __builtin_assume_aligned( q.colptr(l), 32);
    double *Vqj = (double *) __builtin_assume_aligned(Vq.colptr(j), 32);
    double * qj = (double *) __builtin_assume_aligned( q.colptr(j), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        double _sqr_qmi = qm[i]*qm[i];
        ret += sqrt(0.6666666666666666)*(qm[i]*(ql[i]*qm[i]*Vqj[i] + qj[i]*qm[i]*Vql[i] + qj[i]*ql[i]*Vqm[i])*(-3. + 2.*_sqr_qmi) + 3.*V[i]*qj[i]*ql[i]*(-1. + 2.*_sqr_qmi));
    }
    return ret / 3.;
}
double m1_n1_f_m1_l1(int m, int n, int l)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double *Vqn = (double *) __builtin_assume_aligned(Vq.colptr(n), 32);
    double * qn = (double *) __builtin_assume_aligned( q.colptr(n), 32);
    double *Vql = (double *) __builtin_assume_aligned(Vq.colptr(l), 32);
    double * ql = (double *) __builtin_assume_aligned( q.colptr(l), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        ret += 2.*(V[i]*ql[i]*qn[i] + qm[i]*(qm[i]*qn[i]*Vql[i] + ql[i]*qn[i]*Vqm[i] + ql[i]*qm[i]*Vqn[i]));
    }
    return ret / 3.;
}
double m1_n1_f_m1_l2(int m, int n, int l)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double *Vqn = (double *) __builtin_assume_aligned(Vq.colptr(n), 32);
    double * qn = (double *) __builtin_assume_aligned( q.colptr(n), 32);
    double *Vql = (double *) __builtin_assume_aligned(Vq.colptr(l), 32);
    double * ql = (double *) __builtin_assume_aligned( q.colptr(l), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        double _sqr_qli = ql[i]*ql[i];
        ret += V[i]*qn[i]*(-1. + 2.*_sqr_qli) + qm[i]*(2.*ql[i]*qm[i]*qn[i]*Vql[i] + (qn[i]*Vqm[i] + qm[i]*Vqn[i])*(-1. + 2.*_sqr_qli));
    }
    return ret / 3.;
}
double m1_n1_f_m2_l1(int m, int n, int l)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double *Vqn = (double *) __builtin_assume_aligned(Vq.colptr(n), 32);
    double * qn = (double *) __builtin_assume_aligned( q.colptr(n), 32);
    double *Vql = (double *) __builtin_assume_aligned(Vq.colptr(l), 32);
    double * ql = (double *) __builtin_assume_aligned( q.colptr(l), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        double _sqr_qmi = qm[i]*qm[i];
        ret += 4.*V[i]*ql[i]*qm[i]*qn[i] + (qm[i]*qn[i]*Vql[i] + ql[i]*qn[i]*Vqm[i] + ql[i]*qm[i]*Vqn[i])*(-1. + 2.*_sqr_qmi);
    }
    return ret / 3.;
}
double m1_n1_f_m1_l1_j1(int m, int n, int l, int j)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double *Vqn = (double *) __builtin_assume_aligned(Vq.colptr(n), 32);
    double * qn = (double *) __builtin_assume_aligned( q.colptr(n), 32);
    double *Vql = (double *) __builtin_assume_aligned(Vq.colptr(l), 32);
    double * ql = (double *) __builtin_assume_aligned( q.colptr(l), 32);
    double *Vqj = (double *) __builtin_assume_aligned(Vq.colptr(j), 32);
    double * qj = (double *) __builtin_assume_aligned( q.colptr(j), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        ret += V[i]*qj[i]*ql[i]*qn[i] + qm[i]*(ql[i]*qm[i]*qn[i]*Vqj[i] + qj[i]*(qm[i]*qn[i]*Vql[i] + ql[i]*qn[i]*Vqm[i] + ql[i]*qm[i]*Vqn[i]));
    }
    return ret / sqrt(2.);
}
double m1_n2_f_m1_l2(int m, int n, int l)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double *Vqn = (double *) __builtin_assume_aligned(Vq.colptr(n), 32);
    double * qn = (double *) __builtin_assume_aligned( q.colptr(n), 32);
    double *Vql = (double *) __builtin_assume_aligned(Vq.colptr(l), 32);
    double * ql = (double *) __builtin_assume_aligned( q.colptr(l), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        double _sqr_qmi = qm[i]*qm[i];
        double _sqr_qli = ql[i]*ql[i];
        double _sqr_qni = qn[i]*qn[i];
        ret += qn[i]*Vqn[i]*(-1. + 2.*_sqr_qli)*_sqr_qmi + ((V[i] + qm[i]*Vqm[i])*(-1. + 2.*_sqr_qli)*(-1. + 2.*_sqr_qni))/2. + ql[i]*Vql[i]*_sqr_qmi*(-1. + 2.*_sqr_qni);
    }
    return ret / 3.;
}
double m1_n2_f_m2_l1(int m, int n, int l)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double *Vqn = (double *) __builtin_assume_aligned(Vq.colptr(n), 32);
    double * qn = (double *) __builtin_assume_aligned( q.colptr(n), 32);
    double *Vql = (double *) __builtin_assume_aligned(Vq.colptr(l), 32);
    double * ql = (double *) __builtin_assume_aligned( q.colptr(l), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        double _sqr_qmi = qm[i]*qm[i];
        double _sqr_qni = qn[i]*qn[i];
        ret += 4.*V[i]*ql[i]*qm[i]*(-1. + 2.*_sqr_qni) + (-1. + 2.*_sqr_qmi)*(2.*ql[i]*qm[i]*qn[i]*Vqn[i] + (qm[i]*Vql[i] + ql[i]*Vqm[i])*(-1. + 2.*_sqr_qni));
    }
    return ret / 6.;
}
double m1_n2_f_m1_l1_j1(int m, int n, int l, int j)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double *Vqn = (double *) __builtin_assume_aligned(Vq.colptr(n), 32);
    double * qn = (double *) __builtin_assume_aligned( q.colptr(n), 32);
    double *Vql = (double *) __builtin_assume_aligned(Vq.colptr(l), 32);
    double * ql = (double *) __builtin_assume_aligned( q.colptr(l), 32);
    double *Vqj = (double *) __builtin_assume_aligned(Vq.colptr(j), 32);
    double * qj = (double *) __builtin_assume_aligned( q.colptr(j), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        double _sqr_qmi = qm[i]*qm[i];
        double _sqr_qni = qn[i]*qn[i];
        ret += 2.*qj[i]*ql[i]*qn[i]*Vqn[i]*_sqr_qmi + (V[i]*qj[i]*ql[i] + qm[i]*(ql[i]*qm[i]*Vqj[i] + qj[i]*qm[i]*Vql[i] + qj[i]*ql[i]*Vqm[i]))*(-1. + 2.*_sqr_qni);
    }
    return ret / (2.*sqrt(2.));
}
double m2_n1_f_m2_l1(int m, int n, int l)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double *Vqn = (double *) __builtin_assume_aligned(Vq.colptr(n), 32);
    double * qn = (double *) __builtin_assume_aligned( q.colptr(n), 32);
    double *Vql = (double *) __builtin_assume_aligned(Vq.colptr(l), 32);
    double * ql = (double *) __builtin_assume_aligned( q.colptr(l), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        double _sqr_qmi = qm[i]*qm[i];
        ret += qn[i]*(Vql[i] + 2.*(2.*qm[i]*(qm[i]*Vql[i] + ql[i]*Vqm[i])*(-1. + _sqr_qmi) + V[i]*ql[i]*(-1. + 6.*_sqr_qmi))) + ql[i]*Vqn[i]*square(1. - 2.*_sqr_qmi);
    }
    return ret / 6.;
}

double m2_n1_f_m1_l1_j1(int m, int n, int l, int j)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double *Vqn = (double *) __builtin_assume_aligned(Vq.colptr(n), 32);
    double * qn = (double *) __builtin_assume_aligned( q.colptr(n), 32);
    double *Vql = (double *) __builtin_assume_aligned(Vq.colptr(l), 32);
    double * ql = (double *) __builtin_assume_aligned( q.colptr(l), 32);
    double *Vqj = (double *) __builtin_assume_aligned(Vq.colptr(j), 32);
    double * qj = (double *) __builtin_assume_aligned( q.colptr(j), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        double _sqr_qmi = qm[i]*qm[i];
        ret += 4.0*V[i]*qj[i]*ql[i]*qm[i]*qn[i] + (ql[i]*qm[i]*qn[i]*Vqj[i] + qj[i]*(qm[i]*qn[i]*Vql[i] + ql[i]*qn[i]*Vqm[i] + ql[i]*qm[i]*Vqn[i]))*(-1.0 + 2.0*_sqr_qmi);
    }
    return ret / (2.*sqrt(2.0));
}
double m1_n1_p1_f_m1_l1_j1(int m, int n, int p, int l, int j)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double *Vqn = (double *) __builtin_assume_aligned(Vq.colptr(n), 32);
    double * qn = (double *) __builtin_assume_aligned( q.colptr(n), 32);
    double *Vqp = (double *) __builtin_assume_aligned(Vq.colptr(p), 32);
    double * qp = (double *) __builtin_assume_aligned( q.colptr(p), 32);
    double *Vql = (double *) __builtin_assume_aligned(Vq.colptr(l), 32);
    double * ql = (double *) __builtin_assume_aligned( q.colptr(l), 32);
    double *Vqj = (double *) __builtin_assume_aligned(Vq.colptr(j), 32);
    double * qj = (double *) __builtin_assume_aligned( q.colptr(j), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        ret += 4.*(V[i]*qj[i]*ql[i]*qn[i]*qp[i] + qm[i]*(ql[i]*qm[i]*qn[i]*qp[i]*Vqj[i] + qj[i]*(qm[i]*qn[i]*qp[i]*Vql[i] + ql[i]*qn[i]*qp[i]*Vqm[i] + ql[i]*qm[i]*qp[i]*Vqn[i] + ql[i]*qm[i]*qn[i]*Vqp[i])));
    }
    return ret / 5.;
}
double m1_n1_f_m1_n1(int m, int n)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double *Vqn = (double *) __builtin_assume_aligned(Vq.colptr(n), 32);
    double * qn = (double *) __builtin_assume_aligned( q.colptr(n), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        double _sqr_qni = qn[i]*qn[i];
        ret += qm[i]*(V[i]*qm[i] + qn[i]*(qn[i]*Vqm[i] + qm[i]*Vqn[i])) + V[i]*_sqr_qni;
    }
    return ret;
}
double m1_n1_f_m1_n2(int m, int n)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double *Vqn = (double *) __builtin_assume_aligned(Vq.colptr(n), 32);
    double * qn = (double *) __builtin_assume_aligned( q.colptr(n), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        double _sqr_qmi = qm[i]*qm[i];
        double _sqr_qni = qn[i]*qn[i];
        ret += qn[i]*(V[i] + qm[i]*Vqm[i])*(-1. + 2.*_sqr_qni) + _sqr_qmi*(4.*V[i]*qn[i] + Vqn[i]*(-1. + 2.*_sqr_qni));
    }
    return ret / 2.;
}
double m1_n1_f_m2_n1(int m, int n)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double *Vqn = (double *) __builtin_assume_aligned(Vq.colptr(n), 32);
    double * qn = (double *) __builtin_assume_aligned( q.colptr(n), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        double _sqr_qmi = qm[i]*qm[i];
        double _sqr_qni = qn[i]*qn[i];
        ret += (V[i]*qm[i] + qn[i]*(qn[i]*Vqm[i] + qm[i]*Vqn[i]))*(-1. + 2.*_sqr_qmi) + 4.*V[i]*qm[i]*_sqr_qni;
    }
    return ret / 2.;
}
double m1_n1_f_m1_n1_j1(int m, int n, int j)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double *Vqn = (double *) __builtin_assume_aligned(Vq.colptr(n), 32);
    double * qn = (double *) __builtin_assume_aligned( q.colptr(n), 32);
    double *Vqj = (double *) __builtin_assume_aligned(Vq.colptr(j), 32);
    double * qj = (double *) __builtin_assume_aligned( q.colptr(j), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        double _sqr_qni = qn[i]*qn[i];
        ret += 2.*sqrt(2.)*(qm[i]*(V[i]*qj[i]*qm[i] + qn[i]*(qm[i]*qn[i]*Vqj[i] + qj[i]*qn[i]*Vqm[i] + qj[i]*qm[i]*Vqn[i])) + V[i]*qj[i]*_sqr_qni);
    }
    return ret / 3.;
}
double m1_n2_f_m1_n2(int m, int n)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double *Vqn = (double *) __builtin_assume_aligned(Vq.colptr(n), 32);
    double * qn = (double *) __builtin_assume_aligned( q.colptr(n), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        double _sqr_qmi = qm[i]*qm[i];
        double _sqr_qni = qn[i]*qn[i];
        ret += _sqr_qmi*(2.*qn[i]*Vqn[i]*(-1. + _sqr_qni) + V[i]*(-1. + 6.*_sqr_qni)) + ((V[i] + qm[i]*Vqm[i])*square(1. - 2.*_sqr_qni))/2.;
    }
    return ret / 2.;
}
double m1_n2_f_m2_n1(int m, int n)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double *Vqn = (double *) __builtin_assume_aligned(Vq.colptr(n), 32);
    double * qn = (double *) __builtin_assume_aligned( q.colptr(n), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        double _sqr_qmi = qm[i]*qm[i];
        double _sqr_qni = qn[i]*qn[i];
        ret += 4.*V[i]*qm[i]*qn[i]*(-1. + 2.*_sqr_qni) + (-1. + 2.*_sqr_qmi)*(4.*V[i]*qm[i]*qn[i] + (qn[i]*Vqm[i] + qm[i]*Vqn[i])*(-1. + 2.*_sqr_qni));
    }
    return ret / 4.;
}
double m1_n2_f_m1_n1_j1(int m, int n, int j)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double *Vqn = (double *) __builtin_assume_aligned(Vq.colptr(n), 32);
    double * qn = (double *) __builtin_assume_aligned( q.colptr(n), 32);
    double *Vqj = (double *) __builtin_assume_aligned(Vq.colptr(j), 32);
    double * qj = (double *) __builtin_assume_aligned( q.colptr(j), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        double _sqr_qni = qn[i]*qn[i];
        ret += sqrt(2.)*(V[i]*qj[i]*qn[i]*(-1. + 2.*_sqr_qni) + qm[i]*(4.*V[i]*qj[i]*qm[i]*qn[i] + (qm[i]*qn[i]*Vqj[i] + qj[i]*qn[i]*Vqm[i] + qj[i]*qm[i]*Vqn[i])*(-1. + 2.*_sqr_qni)));
    }
    return ret / 3.;
}
double m2_n1_f_m2_n1(int m, int n)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double *Vqn = (double *) __builtin_assume_aligned(Vq.colptr(n), 32);
    double * qn = (double *) __builtin_assume_aligned( q.colptr(n), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        double _sqr_qmi = qm[i]*qm[i];
        double _sqr_qni = qn[i]*qn[i];
        ret += ((V[i] + qn[i]*Vqn[i])*square(1. - 2.*_sqr_qmi))/2. + (2.*qm[i]*Vqm[i]*(-1. + _sqr_qmi) + V[i]*(-1. + 6.*_sqr_qmi))*_sqr_qni;
    }
    return ret / 2.;
}
double m2_n1_f_m1_n1_j1(int m, int n, int j)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double *Vqn = (double *) __builtin_assume_aligned(Vq.colptr(n), 32);
    double * qn = (double *) __builtin_assume_aligned( q.colptr(n), 32);
    double *Vqj = (double *) __builtin_assume_aligned(Vq.colptr(j), 32);
    double * qj = (double *) __builtin_assume_aligned( q.colptr(j), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        double _sqr_qmi = qm[i]*qm[i];
        double _sqr_qni = qn[i]*qn[i];
        ret += sqrt(2.)*((V[i]*qj[i]*qm[i] + qn[i]*(qm[i]*qn[i]*Vqj[i] + qj[i]*qn[i]*Vqm[i] + qj[i]*qm[i]*Vqn[i]))*(-1. + 2.*_sqr_qmi) + 4.*V[i]*qj[i]*qm[i]*_sqr_qni);
    }
    return ret / 3.;
}
double m1_n1_p1_f_m1_n1_j1(int m, int n, int p, int j)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double *Vqn = (double *) __builtin_assume_aligned(Vq.colptr(n), 32);
    double * qn = (double *) __builtin_assume_aligned( q.colptr(n), 32);
    double *Vqp = (double *) __builtin_assume_aligned(Vq.colptr(p), 32);
    double * qp = (double *) __builtin_assume_aligned( q.colptr(p), 32);
    double *Vqj = (double *) __builtin_assume_aligned(Vq.colptr(j), 32);
    double * qj = (double *) __builtin_assume_aligned( q.colptr(j), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        double _sqr_qni = qn[i]*qn[i];
        ret += qm[i]*(V[i]*qj[i]*qm[i]*qp[i] + qn[i]*(qm[i]*qn[i]*qp[i]*Vqj[i] + qj[i]*qn[i]*qp[i]*Vqm[i] + qj[i]*qm[i]*qp[i]*Vqn[i] + qj[i]*qm[i]*qn[i]*Vqp[i])) + V[i]*qj[i]*qp[i]*_sqr_qni;
    }
    return ret;
}
double m1_n1_p1_f_m1_n1_p1(int m, int n, int p)
{
    double *Vqm = (double *) __builtin_assume_aligned(Vq.colptr(m), 32);
    double * qm = (double *) __builtin_assume_aligned( q.colptr(m), 32);
    double *Vqn = (double *) __builtin_assume_aligned(Vq.colptr(n), 32);
    double * qn = (double *) __builtin_assume_aligned( q.colptr(n), 32);
    double *Vqp = (double *) __builtin_assume_aligned(Vq.colptr(p), 32);
    double * qp = (double *) __builtin_assume_aligned( q.colptr(p), 32);
    double ret = 0;
    for (int i=0; i < V.n_rows; i++) {
        double _sqr_qpi = qp[i]*qp[i];
        double _sqr_qmi = qm[i]*qm[i];
        double _sqr_qni = qn[i]*qn[i];
        ret += 4.*((V[i] + qp[i]*Vqp[i])*_sqr_qmi*_sqr_qni + (V[i] + qn[i]*Vqn[i])*_sqr_qmi*_sqr_qpi + (V[i] + qm[i]*Vqm[i])*_sqr_qni*_sqr_qpi);
    }
    return ret / 3.;
}
