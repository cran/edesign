      SUBROUTINE SUBDE1 (DET,A,LDA,NA,AS,LDAS,IND,NI,BS,LDBS,IERR)
      IMPLICIT NONE
*
*  -- user defined driver routine --   
*     
*
*     .. Scalar Arguments ..
      INTEGER            NA, NI, IERR, LDA, LDAS, LDBS
      DOUBLE PRECISION   DET
*     ..
*     .. Array Arguments ..
      INTEGER            IND(*)
      DOUBLE PRECISION   A(LDA,*), AS(LDAS,*), BS(LDBS,*)

*     ..
*
*  Purpose
*  =======
*
*  The function SUBDE1 computes the determinant of the actual submatrix 
*  using the subroutine SUBDET.
*
*
*  Arguments
*  =========
*
*  DET      (output) DOUBLE PRECISION.
*           DET=DET1(1)*10**DET1(2)
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
*           AS stored the actual submatrix.
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
*  IERR     (output) INTEGER
*           the error indicator from DPOFA.  this will be
*           zero unless the matrix being factorized is
*           not positive definite.
*
*
*  =====================================================================
*
*     ..
*     .. Parameters ..
*     ..
*     .. Local Scalars ..
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   DET1(2)
*     ..
*     .. External Functions ..
*     ..
*     .. External Subroutines ..
      EXTERNAL           SUBDET
*     ..
*     .. Intrinsic Functions ..
*     ..
*     .. Executable Statements ..
*
*
      CALL SUBDET(A,LDA,NA,AS,LDAS,IND,NI,BS,LDBS,IERR,DET1)
      DET=DET1(1) * 10 ** DET1(2)
      END
