      SUBROUTINE DSCHUR( N, K, A, LDA, B, LDB, KDET, INFO )
      IMPLICIT NONE
*
*  -- modification of the LAPACK routine DGETF2 --
*     
*     
*     
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LDB, N, K
      DOUBLE PRECISION   KDET
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B(LDB,*)
*     ..
*
*  Purpose
*  =======
*
*  DSCHUR computes the Schur complement of a general m-by-n matrix A
*  using Gaussian elimination.
*
*      [ A[K,K]   A[K,N-K]   ]
*  A = [                     ]
*      [ A[N-K,K] A[N-K,N-K] ]
*  
*  
*  B[K,N-K]=A[N-K,N-K]-A[N-K,K] * A[K,K]^-1 * A[K,N-K]
*  
*  DSCHUR also computes the determinat of the matrix A[K,K].
*
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  K       (input) INTEGER
*          The number of rows of the submatrix A[K,K].  0 <= K <= N
*
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the m by n matrix to be factored.
*          On exit, the factors L and U from the factorization
*          A = P*L*U; the unit diagonal elements of L are not stored.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  B       (output) DOUBLE PRECISION array, dimension (LDB,N)
*          On exit, a N-K by N-K matrix that the elememts of the 
*          Schur complement stored.
*
*
*  LDB     (output) INTEGER
*          The leading dimension of the array B.  0 <= LDB = N-K.
*
*
*  KDET    (output) DOUBLE PRECISION value
*          The determinant of the matrix A[K,K].
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -k, the k-th argument had an illegal value
*
*          > 0: if INFO = k, U(k,k) is exactly zero. 
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO= 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J
      DOUBLE PRECISION   TEN, DET(2)
*     ..
*     .. External Functions ..
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGER, DSCAL, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, DABS
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( K.LT.0 .OR. K.GT.N-1) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LDB.LT.MAX( 1, N-K) ) THEN 
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSCHUR', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. K.EQ.0 .OR. K.EQ.N)
     $   RETURN
*
      DO 10 J = 1, K
*
*     Compute elements J+1:M of J-th column.
*         
         IF( A( J, J ).NE.ZERO ) THEN

            CALL DSCAL( N-J, ONE / A( J, J ), A( J+1, J ), 1 )
*
         ELSE IF( INFO.EQ.0 ) THEN
*
            INFO = J

         END IF
*
*     Update trailing submatrix.
*
         CALL DGER( N-J, N-J, -ONE, A( J+1, J ), 1, A( J, J+1 ), LDA,
     $              A( J+1, J+1 ), LDA )

 10   CONTINUE

      DO 20 I = 1, N-K
         DO 30 J = 1, N-K
*     
*     Store elements A( J, J ) of K+1:N rows and columns in matrix B.
*         
            B( I, J )=A( K+I, K+J )

 30      CONTINUE
 20   CONTINUE

*
*     Determines the determinant of A[K,K].
*     
      DET(1) = 1.0D0
      DET(2) = 0.0D0
      TEN = 10.0D0
      DO 80 I = 1, K
         DET(1) = A(I,I)*DET(1)
         IF (DET(1) .EQ. 0.0D0) GO TO 90
 40      IF (DABS(DET(1)) .GE. 1.0D0) GO TO 50
         DET(1) = TEN*DET(1)
         DET(2) = DET(2) - 1.0D0
         GO TO 40
 50      CONTINUE
 60      IF (DABS(DET(1)) .LT. TEN) GO TO 70
         DET(1) = DET(1)/TEN
         DET(2) = DET(2) + 1.0D0
         GO TO 60
 70      CONTINUE
 80   CONTINUE
 90   CONTINUE
 
      KDET=DET(1)*10**DET(2)

      RETURN
*
*     End of DSCHUR
*
      END
