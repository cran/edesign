      SUBROUTINE SUBDIN(A,LDA,NA,AS,LDAS,IND,NI,BS,LDBS,INV,LDINV,IERR,
     .                  DET)
      IMPLICIT NONE
*
*  -- user defined driver routine --   
*     
*
*     .. Scalar Arguments ..
      INTEGER            NA, NI, IERR, LDA, LDAS, LDINV, LDBS
*     ..
*     .. Array Arguments ..
      INTEGER            IND(*)
      DOUBLE PRECISION   A(LDA,*), AS(LDAS,*), BS(LDBS,*), INV(LDINV,*), 
     .                   DET(*)
*     ..
*
*  Purpose
*  =======
*
*  SUBDIN computes the determinant and the inverse of the actual 
*  submatrix using the Cholesky factorization.
*
*
*  Arguments
*  =========
*
*  A        (input) DOUBLE PRECISION array of DIMENSION (NA,NA).
*           Before entry, the array A must contain matrix of 
*           covariances. 
*           Unchanged on exit.
*
*  NA       (input) INTEGER.
*           On entry, NA specifies the number of rows/columns of A 
*           as declared in the calling (sub) program. NA >= 1.
*           Unchanged on exit.
*
*  AS       (work array) DOUBLE PRECISION, DIMENSION (NA,NA).
*            AS stored the actual submatrix.
*           Changed on exit.
*
*  IND      (input) INTEGER vector, DIMENSION (NA).
*           The first NI elements of the vector correspond to the indices
*           which form the submatrix AS.
*           Unchanged on exit.
*
*  NI       (input) INTEGER.
*           On entry, NI specifies the indices contained in the 
*           actual solution. 0 <= NI <= NA. 
*           Unchanged on exit.
*
*  BS       (work array) DOUBLE PRECISION, DIMENSION (NA,NA).
*           Cholesky decomposition of the matrix A.
*
*  INV      (output) DOUBLE PRECISION array of DIMENSION (NI,NI).
*           On exit, the array INV contains the inverse of the submatrix AS.
*
*  IERR     (output) INTEGER
*           the error indicator from DPOFA.  this will be
*           zero unless the matrix being factorized is
*           not positive definite.
*
*  DET      (output) DOUBLE PRECISION vector, DIMENSION (2).
*           The first argument stores the mantissa and the second
*           gives the exponent of the determinant of the upper 
*           triangular matrix BS.
*
*
*  =====================================================================
*
*     ..
*     .. Parameters ..
*     ..
*     .. Local Scalars ..
      INTEGER            I, J
*     ..
*     .. Local Arrays ..
*     ..
*     .. External Functions ..
*     ..
*     .. External Subroutines ..
      EXTERNAL           PSUBM, CHOL1, CH2INV1, CH2DET
*     ..
*     .. Intrinsic Functions ..
*     ..
*     .. Executable Statements ..
*
*
*     Computing the actual submatrix

      CALL PSUBM(A,LDA,NA,AS,LDAS,IND,NI)

*     Cholesky factorization

      CALL CHOL1(AS,LDA,NI,BS,LDBS,IERR)

*     Determining the inverse  of the submatrix
*     using the output of the Cholesky factorization

      CALL CH2INV1(BS,LDBS,NI,INV,LDINV,IERR)

*     Determining the determinant of the submatrix
*     using the output of the Cholesky factorization

      CALL CH2DET(BS,LDBS,NI,DET,IERR)

      RETURN
      END

