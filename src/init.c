
/*
 * This file contains declarations to interface to the edisgn C and Fortran 
 * functions.
 * 
 */

#include <R.h>
#include <Rinternals.h>

#include "edesign.h"

#include <R_ext/Rdynload.h>

/* Fortran interface descriptions: */

static R_NativePrimitiveArgType grd_t[18] = {
    REALSXP, /* A */
    INTSXP,  /* LDA */
    INTSXP,  /* NA */
    INTSXP,  /* NF */
    INTSXP,  /* NE */
    INTSXP,  /* NS */
    INTSXP,  /* S */
    REALSXP, /* OPT */
    INTSXP,  /* IND */
    REALSXP, /* AS */
    INTSXP,  /* LDAS */
    REALSXP, /* BS */
    INTSXP,  /* LDBS */
    REALSXP, /* INV */
    INTSXP,  /* LDINV */
    REALSXP, /* V */
    REALSXP, /* W */
    INTSXP   /* IERR */
};

static R_NativePrimitiveArgType grddl_t[14] = {
    REALSXP, /* A */
    INTSXP,  /* LDA */
    INTSXP,  /* NA */
    INTSXP,  /* NF */
    INTSXP,  /* NE */
    INTSXP,  /* NS */
    INTSXP,  /* S */
    REALSXP, /* OPT */
    INTSXP,  /* IND */
    REALSXP, /* AS */
    INTSXP,  /* LDAS */
    REALSXP, /* BS */
    INTSXP,  /* LDBS */
    INTSXP   /* IERR */
};

static R_NativePrimitiveArgType change_t[18] = {
    REALSXP, /* A */
    INTSXP,  /* LDA */
    INTSXP,  /* NA */
    INTSXP,  /* NF */
    INTSXP,  /* NE */
    INTSXP,  /* NS */
    INTSXP,  /* S */
    REALSXP, /* OPT */
    INTSXP,  /* IND */
    REALSXP, /* AS */
    INTSXP,  /* LDAS */
    REALSXP, /* BS */
    INTSXP,  /* LDBS */
    REALSXP, /* INV */
    INTSXP,  /* LDINV */
    REALSXP, /* V */
    REALSXP, /* W */
    INTSXP   /* IERR */
};

static R_NativePrimitiveArgType subde1_t[11] = {
    REALSXP, /* DET */
    REALSXP, /* A */
    INTSXP,  /* LDA */
    INTSXP,  /* NA */
    REALSXP, /* AS */
    INTSXP,  /* LDAS */
    INTSXP,  /* IND */
    INTSXP,  /* NI */
    REALSXP, /* BS */
    INTSXP,  /* LDBS */
    INTSXP   /* IERR */
};

static R_FortranMethodDef fortranMethods[] = {
    {"grd",    (DL_FUNC) &F77_SUB(grd), 18, grd_t},       /* greedy */
    {"grddl",  (DL_FUNC) &F77_SUB(grddl), 14, grddl_t},   /* dual.greedy */
    {"change", (DL_FUNC) &F77_SUB(change), 18, change_t}, /* interchange */
    {"subde1", (DL_FUNC) &F77_SUB(subde1), 11, subde1_t}, /* maxentropy */
    {NULL, NULL, 0}
};

/* C interface descriptions: */

static R_NativePrimitiveArgType entrp_t[29] = {
    REALSXP, /* A */
    INTSXP, /* lda */
    INTSXP, /* na */
    INTSXP, /* nf */
    INTSXP, /* ne */
    INTSXP, /* ns */
    INTSXP, /* S */
    REALSXP, /* opt */
    INTSXP, /* S_Work */
    INTSXP, /* F */ 
    INTSXP, /* E */
    INTSXP, /* ind */
    INTSXP, /* ind1 */
    REALSXP, /* As */
    INTSXP, /* ldas */
    REALSXP, /* Bs */
    INTSXP, /* ldbs */
    REALSXP, /* Cs */
    INTSXP, /* ldcs */
    REALSXP, /* Inv */
    INTSXP, /* ldinv */
    REALSXP, /* W */
    REALSXP, /* WORK */
    INTSXP, /* LWORK */
    INTSXP, /* IWORK */
    REALSXP, /* tol */
    INTSXP, /* maxcount */
    INTSXP, /* iter */
    INTSXP /* verbose */
};

static R_CMethodDef cMethods[] = {
    {"entrp", (DL_FUNC) entrp, 29, entrp_t},    /* maxentropy */
    {NULL, NULL, 0}
};


void
R_init_edesign(DllInfo *info)
{
  R_registerRoutines(info, 
		     cMethods, NULL /*callMethods*/, 
		     fortranMethods, NULL/*externalMethods*/);
}
