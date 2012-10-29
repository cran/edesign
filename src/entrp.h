#ifndef ENTRP_H
#define ENTRP_H

#include "R.h" /* for F77_NAME */

/* exported functions, called from R code */

extern void F77_NAME(grd) (double *a,int *lda,int *na,int *nf,int *ne,
                           int *ns,int *s,double *opt,int *ind,double *as,
                           int *ldas,double *bs,int *ldbs,double *inv,
                           int *ldinv,double *v,double *w,int *ierr);
                          
extern void F77_NAME(grddl) (double *a,int *lda,int *na,int *nf,int *ne,
                             int *ns,int *s,double *opt,int *ind,double *as,
                             int *ldas,double *bs,int *ldbs,int *ierr);

extern void F77_NAME(change) (double *a,int *lda,int *na,int *nf,int *ne,
                              int *ns,int *s,double *opt,int *ind,double *as,
                              int *ldas,double *bs,int *ldbs,double *inv,
                              int *ldinv,double *v,double *w,int *ierr);

extern void F77_NAME(subde1) (double *det,double *a,int *lda,int *na,
                              double *as,int *ldas,int *ind,int *ni,
                              double *bs,int *ldbs,int *ierr);

/* internal functions, called from C code */

extern double F77_NAME(upbnd)(double*,int*,int*,int*,int*,int*,int*,double*,
                              int*,double*,int*,double*,int*,double*,int*,
                              int*,int*,double*,double*,int*,int*);

extern double F77_NAME(subdet)(double*,int*,int*,double*,int*,int*, double*,
                               int*,int*,double*);

extern double F77_NAME(chdet)(double*,int*,int*);

extern void F77_NAME(chol1)(double*,int*,int*,double*,int*,int*);

extern void F77_NAME(ivecpr)(int*,int*); 

extern void F77_NAME(psubm)(double*,int*,int*,double*,int*,int*,int*); 


/* main C function */

void entrp(double *A,int *lda,int *na,int *nf,int *ne,int *ns,int *S,
           double *opt,int *S_Work,int *F,int *E,int *ind,int *ind1,
           double *As,int *ldas,double *Bs,int *ldbs,double *Cs,int *ldcs,
           double *Inv,int *ldinv,double *W,double* WORK,int *LWORK,
           int *IWORK,double *tol,int *maxcount,int *iter,int *verbose);


#endif
