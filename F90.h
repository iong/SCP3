/*
 *  F90.h
 *  SCP3
 *
 *  Created by Ionuț Georgescu on 4/2/13.
 *  Copyright 2013 UCI Chemistry. All rights reserved.
 *
 */


extern "C" {
    void sobol_stdnormal_c(int64_t d, int64_t *skip, void *x);
    void TIP4P_UF(int N, double *r, double *U, double *UX);
    void dsyevr_(char *, char *, char *,  int *,  double *,  int *,  double *,
                 double* ,  int *IL,  int *IU, double *ABSTOL,  int *M,
                 double *W, double *Z, int *LDZ, int *ISUPPZ, double *WORK,
                 int *LWORK, int *IWORK, int *LIWORK, int *INFO);
    double dlamch_(char *);

    void whbb_pes_init(int nw);
    void whbb_fgrad(int nw, double *x, double eps, double *V, double *Vx);
}
