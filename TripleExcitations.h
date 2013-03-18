double f_k1(int k)
{
    return Vq[k]/sqrt(2.);
}
double f_k2(int k)
{
    return (q[k]*Vq[k])/sqrt(2.);
}
double f_k3(int k)
{
    return (4.*V*q[k] + Vq[k]*(-3. + 2.*sqr(q[k])))/(2.*sqrt(3.));
}
double f_k1_l1(int k, int l)
{
    return (q[l]*Vq[k] + q[k]*Vq[l])/2.;
}
double f_k1_l2(int k, int l)
{
    return (q[k]*q[l]*Vq[l] + Vq[k]*(-0.5 + sqr(q[l])))/2.;
}
double f_k2_l1(int k, int l)
{
    return (q[k]*q[l]*Vq[k] + Vq[l]*(-0.5 + sqr(q[k])))/2.;
}
double f_k1_l1_j1(int k, int l, int j)
{
    return (sqrt(2.)*(q[k]*q[l]*Vq[j] + q[j]*q[l]*Vq[k] + q[j]*q[k]*Vq[l]))/3.;
}

double m1_f_k1(int m, int k)
{
    return (q[m]*Vq[k] + q[k]*Vq[m])/2.;
}

double m1_f_k2(int m, int k)
{
    return (q[k]*q[m]*Vq[k] + Vq[m]*(-0.5 + sqr(q[k])))/2.;
}


double m1_f_k3(int m, int k)
{
    return (4.*V*q[k]*q[m] + (q[m]*Vq[k] + q[k]*Vq[m])*(-3. + 2.*sqr(q[k])))/(2.*sqrt(6.));
}


double m1_f_k1_l1(int m, int k, int l)
{
    return (sqrt(2.)*(q[l]*q[m]*Vq[k] + q[k]*q[m]*Vq[l] + q[k]*q[l]*Vq[m]))/3.;
}


double m1_f_k1_l2(int m, int k, int l)
{
    return (2.*q[k]*q[l]*q[m]*Vq[l] + q[m]*Vq[k]*(-1. + 2.*sqr(q[l])) + q[k]*Vq[m]*(-1. + 2.*sqr(q[l])))/(3.*sqrt(2.));
}
double m1_f_k2_l1(int m, int k, int l)
{
    return (2.*q[k]*q[l]*q[m]*Vq[k] + (q[m]*Vq[l] + q[l]*Vq[m])*(-1. + 2.*sqr(q[k])))/(3.*sqrt(2.));
}
double m1_f_k1_l1_j1(int m, int k, int l, int j)
{
    return (q[k]*q[l]*q[m]*Vq[j] + q[j]*(q[l]*q[m]*Vq[k] + q[k]*q[m]*Vq[l] + q[k]*q[l]*Vq[m]))/2.;
}
double m2_f_k2(int m, int k)
{
    return (q[m]*Vq[m]*(-1. + 2.*sqr(q[k])) + q[k]*Vq[k]*(-1. + 2.*sqr(q[m])))/4.;
}
double m2_f_k3(int m, int k)
{
    return (4.*V*q[k]*(-1. + 2.*sqr(q[m])) + (-3. + 2.*sqr(q[k]))*(2.*q[k]*q[m]*Vq[m] + Vq[k]*(-1. + 2.*sqr(q[m]))))/(4.*sqrt(6.));
}
double m2_f_k1_l1(int m, int k, int l)
{
    return (2.*q[k]*q[l]*q[m]*Vq[m] + (q[l]*Vq[k] + q[k]*Vq[l])*(-1. + 2.*sqr(q[m])))/(3.*sqrt(2.));
}
double m2_f_k1_l2(int m, int k, int l)
{
    return (2.*q[k]*q[m]*Vq[m]*(-1. + 2.*sqr(q[l])) + (2.*q[k]*q[l]*Vq[l] + Vq[k]*(-1. + 2.*sqr(q[l])))*(-1. + 2.*sqr(q[m])))/(6.*sqrt(2.));
}
double m2_f_k2_l1(int m, int k, int l)
{
    return (2.*q[k]*q[l]*Vq[k]*(-1. + 2.*sqr(q[m])) + (-1. + 2.*sqr(q[k]))*(2.*q[l]*q[m]*Vq[m] + Vq[l]*(-1. + 2.*sqr(q[m]))))/(6.*sqrt(2.));
}
double m2_f_k1_l1_j1(int m, int k, int l, int j)
{
    return (2.*q[j]*q[k]*q[l]*q[m]*Vq[m] + (q[k]*q[l]*Vq[j] + q[j]*q[l]*Vq[k] + q[j]*q[k]*Vq[l])*(-1. + 2.*sqr(q[m])))/4.;
}
double m3_f_k3(int m, int k)
{
    return (4.*V*q[k]*q[m]*(-3. + 2.*sqr(q[m])) + (-3. + 2.*sqr(q[k]))*(4.*V*q[k]*q[m] + (q[m]*Vq[k] + q[k]*Vq[m])*(-3. + 2.*sqr(q[m]))))/12.;
}
double m3_f_k1_l1(int m, int k, int l)
{
    return (4.*V*q[k]*q[l]*q[m] + (q[l]*q[m]*Vq[k] + q[k]*q[m]*Vq[l] + q[k]*q[l]*Vq[m])*(-3. + 2.*sqr(q[m])))/(3.*sqrt(3.));
}
double m3_f_k1_l2(int m, int k, int l)
{
    return (4.*V*q[k]*q[m]*(-1. + 2.*sqr(q[l])) + (2.*q[k]*q[l]*q[m]*Vq[l] + q[m]*Vq[k]*(-1. + 2.*sqr(q[l])) + q[k]*Vq[m]*(-1. + 2.*sqr(q[l])))*(-3. + 2.*sqr(q[m])))/(6.*sqrt(3.));
}
double m3_f_k2_l1(int m, int k, int l)
{
    return (4.*V*q[l]*q[m]*(-1. + 2.*sqr(q[k])) + (2.*q[k]*q[l]*q[m]*Vq[k] + (q[m]*Vq[l] + q[l]*Vq[m])*(-1. + 2.*sqr(q[k])))*(-3. + 2.*sqr(q[m])))/(6.*sqrt(3.));
}
double m3_f_k1_l1_j1(int m, int k, int l, int j)
{
    return (4.*V*q[j]*q[k]*q[l]*q[m] + (q[k]*q[l]*q[m]*Vq[j] + q[j]*(q[l]*q[m]*Vq[k] + q[k]*q[m]*Vq[l] + q[k]*q[l]*Vq[m]))*(-3. + 2.*sqr(q[m])))/(2.*sqrt(6.));
}
double m1_n1_f_k1_l1(int m, int n, int k, int l)
{
    return (q[l]*q[m]*q[n]*Vq[k] + q[k]*(q[m]*q[n]*Vq[l] + q[l]*q[n]*Vq[m] + q[l]*q[m]*Vq[n]))/2.;
}
double m1_n1_f_k1_l2(int m, int n, int k, int l)
{
    return (2.*q[k]*q[l]*q[m]*q[n]*Vq[l] + q[m]*q[n]*Vq[k]*(-1. + 2.*sqr(q[l])) + q[k]*(q[n]*Vq[m] + q[m]*Vq[n])*(-1. + 2.*sqr(q[l])))/4.;
}
double m1_n1_f_k2_l1(int m, int n, int k, int l)
{
    return (2.*q[k]*q[l]*q[m]*q[n]*Vq[k] + (q[m]*q[n]*Vq[l] + q[l]*q[n]*Vq[m] + q[l]*q[m]*Vq[n])*(-1. + 2.*sqr(q[k])))/4.;
}
double m1_n1_f_k1_l1_j1(int m, int n, int k, int l, int j)
{
    return (2.*sqrt(2.)*(q[k]*q[l]*q[m]*q[n]*Vq[j] + q[j]*(q[l]*q[m]*q[n]*Vq[k] + q[k]*(q[m]*q[n]*Vq[l] + q[l]*q[n]*Vq[m] + q[l]*q[m]*Vq[n]))))/5.;
}
double m1_n2_f_k1_l2(int m, int n, int k, int l)
{
    return (2.*q[k]*q[m]*q[n]*Vq[n]*(-1. + 2.*sqr(q[l])) + (2.*q[k]*q[l]*q[m]*Vq[l] + q[m]*Vq[k]*(-1. + 2.*sqr(q[l])) + q[k]*Vq[m]*(-1. + 2.*sqr(q[l])))*(-1. + 2.*sqr(q[n])))/8.;
}
double m1_n2_f_k2_l1(int m, int n, int k, int l)
{
    return (2.*q[k]*q[l]*q[m]*Vq[k]*(-1. + 2.*sqr(q[n])) + (-1. + 2.*sqr(q[k]))*(2.*q[l]*q[m]*q[n]*Vq[n] + (q[m]*Vq[l] + q[l]*Vq[m])*(-1. + 2.*sqr(q[n]))))/8.;
}
double m1_n2_f_k1_l1_j1(int m, int n, int k, int l, int j)
{
    return (sqrt(2.)*(2.*q[j]*q[k]*q[l]*q[m]*q[n]*Vq[n] + (q[k]*q[l]*q[m]*Vq[j] + q[j]*(q[l]*q[m]*Vq[k] + q[k]*q[m]*Vq[l] + q[k]*q[l]*Vq[m]))*(-1. + 2.*sqr(q[n]))))/5.;
}
double m2_n1_f_k2_l1(int m, int n, int k, int l)
{
    return (2.*q[k]*q[l]*q[n]*Vq[k]*(-1. + 2.*sqr(q[m])) + (-1. + 2.*sqr(q[k]))*(2.*q[l]*q[m]*q[n]*Vq[m] + q[n]*Vq[l]*(-1. + 2.*sqr(q[m])) + q[l]*Vq[n]*(-1. + 2.*sqr(q[m]))))/8.;
}
double m2_n1_f_k1_l1_j1(int m, int n, int k, int l, int j)
{
    return (sqrt(2.)*(2.*q[j]*q[k]*q[l]*q[m]*q[n]*Vq[m] + q[n]*(q[k]*q[l]*Vq[j] + q[j]*q[l]*Vq[k] + q[j]*q[k]*Vq[l])*(-1. + 2.*sqr(q[m])) + q[j]*q[k]*q[l]*Vq[n]*(-1. + 2.*sqr(q[m]))))/5.;
}
double m1_n1_p1_f_k1_l1_j1(int m, int n, int p, int k, int l, int j)
{
    return (2.*(q[k]*q[l]*q[m]*q[n]*q[p]*Vq[j] + q[j]*q[l]*q[m]*q[n]*q[p]*Vq[k] + q[j]*q[k]*(q[m]*q[n]*q[p]*Vq[l] + q[l]*q[n]*q[p]*Vq[m] + q[l]*q[m]*q[p]*Vq[n] + q[l]*q[m]*q[n]*Vq[p])))/3.;
}
double m1_f_m1(int m)
{
    return V + q[m]*Vq[m];
}
double m1_f_m2(int m)
{
    return 2.*V*q[m] + Vq[m]*(-0.5 + sqr(q[m]));
}
double m1_f_m3(int m)
{
    return (q[m]*Vq[m]*(-3. + 2.*sqr(q[m])) + V*(-3. + 6.*sqr(q[m])))/sqrt(6.);
}
double m1_f_m1_l1(int m, int l)
{
    return (V*q[l] + q[m]*(q[m]*Vq[l] + q[l]*Vq[m]))/sqrt(2.);
}
double m1_f_m1_l2(int m, int l)
{
    return ((V + q[m]*Vq[m])*(-1. + 2.*sqr(q[l])) + 2.*q[l]*Vq[l]*sqr(q[m]))/(2.*sqrt(2.));
}
double m1_f_m2_l1(int m, int l)
{
    return (4.*V*q[l]*q[m] + (q[m]*Vq[l] + q[l]*Vq[m])*(-1. + 2.*sqr(q[m])))/(2.*sqrt(2.));
}
double m1_f_m1_l1_j1(int m, int l, int j)
{
    return (2.*(V*q[j]*q[l] + q[m]*(q[l]*q[m]*Vq[j] + q[j]*q[m]*Vq[l] + q[j]*q[l]*Vq[m])))/3.;
}
double m2_f_m2(int m)
{
    return q[m]*Vq[m]*(-1. + sqr(q[m])) + V*(-0.5 + 3.*sqr(q[m]));
}

double m2_f_m3(int m)
{
    return (Vq[m]*(3. + 4.*sqr(sqr(q[m])) - 8.*sqr(q[m])) + 16.*V*q[m]*(-1. + sqr(q[m])))/(2.*sqrt(6.));
}
double m2_f_m1_l1(int m, int l)
{
    return (4.*V*q[l]*q[m] + (q[m]*Vq[l] + q[l]*Vq[m])*(-1. + 2.*sqr(q[m])))/(2.*sqrt(2.));
}
double m2_f_m1_l2(int m, int l)
{
    return (4.*V*q[m]*(-1. + 2.*sqr(q[l])) + (2.*q[l]*q[m]*Vq[l] + Vq[m]*(-1. + 2.*sqr(q[l])))*(-1. + 2.*sqr(q[m])))/(4.*sqrt(2.));
}
double m2_f_m2_l1(int m, int l)
{
    return (Vq[l] + 2.*(2.*q[m]*(q[m]*Vq[l] + q[l]*Vq[m])*(-1. + sqr(q[m])) + V*q[l]*(-1. + 6.*sqr(q[m]))))/(4.*sqrt(2.));
}
double m2_f_m1_l1_j1(int m, int l, int j)
{
    return (4.*V*q[j]*q[l]*q[m] + (q[l]*q[m]*Vq[j] + q[j]*q[m]*Vq[l] + q[j]*q[l]*Vq[m])*(-1. + 2.*sqr(q[m])))/3.;
}
double m3_f_m3(int m)
{
    return ((-3. + 2.*sqr(q[m]))*(q[m]*Vq[m]*(-3. + 2.*sqr(q[m])) + V*(-3. + 1.0*sqr(q[m]))))/6.;
}
double m3_f_m1_l1(int m, int l)
{
    return (q[m]*(q[m]*Vq[l] + q[l]*Vq[m])*(-3. + 2.*sqr(q[m])) + 3.*V*q[l]*(-1. + 2.*sqr(q[m])))/(2.*sqrt(3.));
}
double m3_f_m1_l2(int m, int l)
{
    return (q[m]*(2.*q[l]*q[m]*Vq[l] + Vq[m]*(-1. + 2.*sqr(q[l])))*(-3. + 2.*sqr(q[m])) + 3.*V*(-1. + 2.*sqr(q[l]))*(-1. + 2.*sqr(q[m])))/(4.*sqrt(3.));
}
double m3_f_m2_l1(int m, int l)
{
    return ((q[m]*Vq[l] + q[l]*Vq[m])*(3. + 4.*sqr(sqr(q[m])) - 8.*sqr(q[m])) + 16.*V*q[l]*q[m]*(-1. + sqr(q[m])))/(4.*sqrt(3.));
}
double m3_f_m1_l1_j1(int m, int l, int j)
{
    return (sqrt(0.6666666666666666)*(q[m]*(q[l]*q[m]*Vq[j] + q[j]*q[m]*Vq[l] + q[j]*q[l]*Vq[m])*(-3. + 2.*sqr(q[m])) + 3.*V*q[j]*q[l]*(-1. + 2.*sqr(q[m]))))/3.;
}
double m1_n1_f_m1_l1(int m, int n, int l)
{
    return (2.*(V*q[l]*q[n] + q[m]*(q[m]*q[n]*Vq[l] + q[l]*q[n]*Vq[m] + q[l]*q[m]*Vq[n])))/3.;
}
double m1_n1_f_m1_l2(int m, int n, int l)
{
    return (V*q[n]*(-1. + 2.*sqr(q[l])) + q[m]*(2.*q[l]*q[m]*q[n]*Vq[l] + (q[n]*Vq[m] + q[m]*Vq[n])*(-1. + 2.*sqr(q[l]))))/3.;
}
double m1_n1_f_m2_l1(int m, int n, int l)
{
    return (4.*V*q[l]*q[m]*q[n] + (q[m]*q[n]*Vq[l] + q[l]*q[n]*Vq[m] + q[l]*q[m]*Vq[n])*(-1. + 2.*sqr(q[m])))/3.;
}
double m1_n1_f_m1_l1_j1(int m, int n, int l, int j)
{
    return (V*q[j]*q[l]*q[n] + q[m]*(q[l]*q[m]*q[n]*Vq[j] + q[j]*(q[m]*q[n]*Vq[l] + q[l]*q[n]*Vq[m] + q[l]*q[m]*Vq[n])))/sqrt(2.);
}
double m1_n2_f_m1_l2(int m, int n, int l)
{
    return (q[n]*Vq[n]*(-1. + 2.*sqr(q[l]))*sqr(q[m]) + ((V + q[m]*Vq[m])*(-1. + 2.*sqr(q[l]))*(-1. + 2.*sqr(q[n])))/2. + q[l]*Vq[l]*sqr(q[m])*(-1. + 2.*sqr(q[n])))/3.;
}
double m1_n2_f_m2_l1(int m, int n, int l)
{
    return (4.*V*q[l]*q[m]*(-1. + 2.*sqr(q[n])) + (-1. + 2.*sqr(q[m]))*(2.*q[l]*q[m]*q[n]*Vq[n] + (q[m]*Vq[l] + q[l]*Vq[m])*(-1. + 2.*sqr(q[n]))))/6.;
}
double m1_n2_f_m1_l1_j1(int m, int n, int l, int j)
{
    return (2.*q[j]*q[l]*q[n]*Vq[n]*sqr(q[m]) + (V*q[j]*q[l] + q[m]*(q[l]*q[m]*Vq[j] + q[j]*q[m]*Vq[l] + q[j]*q[l]*Vq[m]))*(-1. + 2.*sqr(q[n])))/(2.*sqrt(2.));
}
double m2_n1_f_m2_l1(int m, int n, int l)
{
    return (q[n]*(Vq[l] + 2.*(2.*q[m]*(q[m]*Vq[l] + q[l]*Vq[m])*(-1. + sqr(q[m])) + V*q[l]*(-1. + 6.*sqr(q[m])))) + q[l]*Vq[n]*sqr(1. - 2.*sqr(q[m])))/6.;
}

#include <cstdio>

double m2_n1_f_m1_l1_j1(int m, int n, int l, int j)
{
    return (4.0*V*q[j]*q[l]*q[m]*q[n] + (q[l]*q[m]*q[n]*Vq[j] + q[j]*(q[m]*q[n]*Vq[l] + q[l]*q[n]*Vq[m] + q[l]*q[m]*Vq[n]))*(-1.0 + 2.0*sqr(q[m])))/(2.*sqrt(2.0));
}
double m1_n1_p1_f_m1_l1_j1(int m, int n, int p, int l, int j)
{
    return (4.*(V*q[j]*q[l]*q[n]*q[p] + q[m]*(q[l]*q[m]*q[n]*q[p]*Vq[j] + q[j]*(q[m]*q[n]*q[p]*Vq[l] + q[l]*q[n]*q[p]*Vq[m] + q[l]*q[m]*q[p]*Vq[n] + q[l]*q[m]*q[n]*Vq[p]))))/5.;
}
double m1_n1_f_m1_n1(int m, int n)
{
    return q[m]*(V*q[m] + q[n]*(q[n]*Vq[m] + q[m]*Vq[n])) + V*sqr(q[n]);
}
double m1_n1_f_m1_n2(int m, int n)
{
    return (q[n]*(V + q[m]*Vq[m])*(-1. + 2.*sqr(q[n])) + sqr(q[m])*(4.*V*q[n] + Vq[n]*(-1. + 2.*sqr(q[n]))))/2.;
}
double m1_n1_f_m2_n1(int m, int n)
{
    return ((V*q[m] + q[n]*(q[n]*Vq[m] + q[m]*Vq[n]))*(-1. + 2.*sqr(q[m])) + 4.*V*q[m]*sqr(q[n]))/2.;
}
double m1_n1_f_m1_n1_j1(int m, int n, int j)
{
    return (2.*sqrt(2.)*(q[m]*(V*q[j]*q[m] + q[n]*(q[m]*q[n]*Vq[j] + q[j]*q[n]*Vq[m] + q[j]*q[m]*Vq[n])) + V*q[j]*sqr(q[n])))/3.;
}
double m1_n2_f_m1_n2(int m, int n)
{
    return (sqr(q[m])*(2.*q[n]*Vq[n]*(-1. + sqr(q[n])) + V*(-1. + 6.*sqr(q[n]))) + ((V + q[m]*Vq[m])*sqr(1. - 2.*sqr(q[n])))/2.)/2.;
}
double m1_n2_f_m2_n1(int m, int n)
{
    return (4.*V*q[m]*q[n]*(-1. + 2.*sqr(q[n])) + (-1. + 2.*sqr(q[m]))*(4.*V*q[m]*q[n] + (q[n]*Vq[m] + q[m]*Vq[n])*(-1. + 2.*sqr(q[n]))))/4.;
}
double m1_n2_f_m1_n1_j1(int m, int n, int j)
{
    return (sqrt(2.)*(V*q[j]*q[n]*(-1. + 2.*sqr(q[n])) + q[m]*(4.*V*q[j]*q[m]*q[n] + (q[m]*q[n]*Vq[j] + q[j]*q[n]*Vq[m] + q[j]*q[m]*Vq[n])*(-1. + 2.*sqr(q[n])))))/3.;
}
double m2_n1_f_m2_n1(int m, int n)
{
    return (((V + q[n]*Vq[n])*sqr(1. - 2.*sqr(q[m])))/2. + (2.*q[m]*Vq[m]*(-1. + sqr(q[m])) + V*(-1. + 6.*sqr(q[m])))*sqr(q[n]))/2.;
}
double m2_n1_f_m1_n1_j1(int m, int n, int j)
{
    return (sqrt(2.)*((V*q[j]*q[m] + q[n]*(q[m]*q[n]*Vq[j] + q[j]*q[n]*Vq[m] + q[j]*q[m]*Vq[n]))*(-1. + 2.*sqr(q[m])) + 4.*V*q[j]*q[m]*sqr(q[n])))/3.;
}
double m1_n1_p1_f_m1_n1_j1(int m, int n, int p, int j)
{
    return q[m]*(V*q[j]*q[m]*q[p] + q[n]*(q[m]*q[n]*q[p]*Vq[j] + q[j]*q[n]*q[p]*Vq[m] + q[j]*q[m]*q[p]*Vq[n] + q[j]*q[m]*q[n]*Vq[p])) + V*q[j]*q[p]*sqr(q[n]);
}
double m1_n1_p1_f_m1_n1_p1(int m, int n, int p)
{
    return (4.*((V + q[p]*Vq[p])*sqr(q[m])*sqr(q[n]) + (V + q[n]*Vq[n])*sqr(q[m])*sqr(q[p]) + (V + q[m]*Vq[m])*sqr(q[n])*sqr(q[p])))/3.;
}
