      SUBROUTINE CH2DET(A,LDA,N,DET,IERR)
      IMPLICIT NONE
*
*  -- user defined driver routine --   
*     
*
*     .. Scalar Arguments ..
      INTEGER            LDA, N, IERR
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A(LDA,*), DET(*)
*     ..
*
*  Purpose
*  =======
*
*  CH2DET computes the determinant of a matrix proceeding from 
*  the Cholesky factorization A of this matrix.
*
*
*  Arguments
*  =========
*
*  A        (input) DOUBLE PRECISION array of DIMENSION (N,N).
*           Before entry, the array A must contain in its upper half 
*           the Cholesky factor of the matrix A.
*           Unchanged on exit.
*
*  LDA      (input) INTEGER.
*           LDA is the leading dimension of the array A.
*           Unchanged on exit.
*
*  N        (input) INTEGER.
*           N is the order of the matrix A.
*           Unchanged on exit.
*
*  DET      (output) DOUBLE PRECISION vector, DIMENSION (2).
*           The first argument stores the mantissa and the second
*           gives the exponent of the determinant of the upper 
*           triangular matrix BS.
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
      INTEGER            I
      DOUBLE PRECISION   TEN
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
      DET(1) = 1.0D0
      DET(2) = 0.0D0
      TEN = 10.0D0
      DO 50 I = 1, N
         DET(1) = A(I,I)*A(I,I)*DET(1)
         IF (DET(1) .EQ. 0.0D0) GO TO 60
 10      IF (DABS(DET(1)) .GE. 1.0D0) GO TO 20
         DET(1) = TEN*DET(1)
         DET(2) = DET(2) - 1.0D0
         GO TO 10
 20      CONTINUE
 30      IF (DABS(DET(1)) .LT. TEN) GO TO 40
         DET(1) = DET(1)/TEN
         DET(2) = DET(2) + 1.0D0
         GO TO 30
 40      CONTINUE
 50   CONTINUE
 60   CONTINUE
      IERR=0
      RETURN
      END

