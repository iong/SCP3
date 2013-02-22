//
//  DoubleExcitations.h
//  SCP_Double_Excitations
//
//  Created by Ionut Georgescu on 2/21/13.
//
//

#ifndef SCP_Double_Excitations_DoubleExcitations_h
#define SCP_Double_Excitations_DoubleExcitations_h

double m0_f_k1(int k)
{
    return Vq[k]/sqrt(2);
}
double m0_f_k2(int k)
{
    return (q[k]*Vq[k])/sqrt(2);
}
double m1_f_k1(int m, int k)
{
    return (q[m]*Vq[k] + q[k]*Vq[m])/2.;
}
double m1_f_k2(int m, int k)
{
    return (q[k]*q[m]*Vq[k] + Vq[m]*(-0.5 + sqr(q[k])))/2.;
}
double m2_f_k2(int m, int k)
{
    return (q[m]*Vq[m]*(-1 + 2*sqr(q[k])) + q[k]*Vq[k]*(-1 + 2*sqr(q[m])))/4.;
}
double m1_f_m1(int m)
{
    return V + q[m]*Vq[m];
}
double m1_f_m2(int m)
{
    return 2*V*q[m] + Vq[m]*(-0.5 + sqr(q[m]));
}
double m2_f_m2(int m)
{
    return q[m]*Vq[m]*(-1 + sqr(q[m])) + V*(-0.5 + 3*sqr(q[m]));
}
double m0_f_k1_l1(int k, int l)
{
    return (2*V*q[k]*q[l] + q[l]*Vq[k] + q[k]*Vq[l])/3.;
}
double m1_f_k1_l1(int m, int k, int l)
{
    return ((q[l]*q[m]*Vq[k] + q[k]*q[m]*Vq[l] + q[k]*q[l]*Vq[m])*sqrt(2.))/3.;
}
double m2_f_k1_l1(int m, int k, int l)
{
    return (2*q[k]*q[l]*q[m]*Vq[m] + (q[l]*Vq[k] + q[k]*Vq[l])*(-1 + 2*sqr(q[m])))/(3.*sqrt(2));
}
double m1_f_m1_l1(int m, int l)
{
    return (V*q[l] + q[m]*(q[m]*Vq[l] + q[l]*Vq[m]))/sqrt(2);
}
double m2_f_m1_l1(int m, int l)
{
    return (4*V*q[l]*q[m] + (q[m]*Vq[l] + q[l]*Vq[m])*(-1 + 2*sqr(q[m])))/(2.*sqrt(2));
}
double m1_n1_f_k1_l1(int m, int n, int k, int l)
{
    return (q[l]*q[m]*q[n]*Vq[k] + q[k]*(q[m]*q[n]*Vq[l] + q[l]*q[n]*Vq[m] + q[l]*q[m]*Vq[n]))/2.;
}
double m1_n1_f_m1_l1(int m, int n, int l)
{
    return (2*(V*q[l]*q[n] + q[m]*(q[m]*q[n]*Vq[l] + q[l]*q[n]*Vq[m] + q[l]*q[m]*Vq[n])))/3.;
}
double m1_n1_f_m1_n1(int m, int n)
{
    return q[m]*(V*q[m] + q[n]*(q[n]*Vq[m] + q[m]*Vq[n])) + V*sqr(q[n]);
}


#endif
