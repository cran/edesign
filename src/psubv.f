      SUBROUTINE PSUBV(A,LDA,NA,V,IND,NI,I)
      IMPLICIT NONE

*
*  -- user defined driver routine --   
*     
*
*     .. Scalar Arguments ..
      INTEGER            NA, NI, I, LDA
*     ..
*     .. Array Arguments ..
      INTEGER            IND(*)
      DOUBLE PRECISION   A(LDA,*), V(*)
*     ..
*
*  Purpose
*  =======
*
*  PSUBV computes the I-th column vector V of A. The elements of V 
*  are reduced to the actual index set IND.
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
*  V        (output) DOUBLE PRECISION vector of DIMENSION (NI).
*           On exit, the I-th column vector of A reduced to the 
*           actual index set IND. 
*
*  IND      (input) INTEGER vector, DIMENSION (NI).
*           On entry, IND specifies the indices forced into 
*           the actual solution.
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
      INTEGER            J
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
      DO 10 J=1,NI
         V(J)=0
         V(J)=A(IND(J),I)
  10   CONTINUE
      END
