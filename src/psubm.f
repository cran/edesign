      SUBROUTINE PSUBM(A,LDA,NA,AS,LDAS,IND,NI)
      IMPLICIT NONE
*
*  -- user defined driver routine --   
*     
*
*     .. Scalar Arguments ..
      INTEGER            NA, NI, LDA, LDAS
*     ..
*     .. Array Arguments ..
      INTEGER            IND(*)
      DOUBLE PRECISION   A(LDA,*), AS(LDAS,*)
*     ..
*
*  Purpose
*  =======
*
*  PSUBM computes proceeding from the actual index set IND the 
*  principal submatrix AS of A.
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
*           as declared in the calling (sub) program. NA >= 0.
*           Unchanged on exit.
*
*  AS       (output) DOUBLE PRECISION array of DIMENSION (NI,NI).
*           On exit, the principal submatrix AS. 
*
*  IND      (input) INTEGER vector, DIMENSION (NI).
*           On entry, IND specifies the indices forced into 
*           the principal submatrix.
*           Unchanged on exit.
*
*  NI       (input) INTEGER.
*           On entry, NI specifies the dimension of IND. 0 <= NI <= NA. 
*           Unchanged on exit.
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
*     ..
*     .. Intrinsic Functions ..
*     ..
*     .. Executable Statements ..
*
*
      DO 10 I=1,NI
         DO 20 J=1,NI
            AS(I,J)=A(IND(I),IND(J))
 20      CONTINUE
 10   CONTINUE
      END
