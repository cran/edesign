      DOUBLE PRECISION FUNCTION UPBND(A,LDA,NA,F,E,NE,NS,AS,LDAS,BS,
     .                                LDBS,CS,LDCS,INV,LDINV,IND,
     .                                IERR,W,WORK,LWORK,IWORK)
      IMPLICIT NONE
*
*  -- user defined driver routine --   
*     
*
*     .. Scalar Arguments ..
      INTEGER            NA, NE, NS, IERR, LWORK, LDA, LDAS, LDBS, LDCS, 
     .                   LDINV
*     ..
*     .. Array Arguments ..
      INTEGER            F(*), E(*), IND(*), IWORK(*)
      DOUBLE PRECISION   A(LDA,*), AS(LDAS,*), CS(LDCS, *), BS(LDBS,*),
     .                   INV(LDINV,*), W(*), WORK(*)
*     ..
*
*  Purpose
*  =======
*
*  UPBND is a function to establish the upper bounds for the entropy 
*  required by the branch and bound strategy described in 
*  Ko, Lee, Queyranne (1994), An exact algorithm for maximum entropy 
*  sampling, Operations Research 43 (1995), 684-691.       
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
*  F        (input) INTEGER vector of DIMENSION (NE).
*           On entry, F specifies the additional indices are forced  
*           into the actual solution. 
*           Unchanged on exit.
*
*  E        (input) INTEGER vector of DIMENSION (NE).
*           On entry, E specifies the indices still treated  
*           as eligible for the actual solution. 
*           Unchanged on exit.
*
*  NE       (input) INTEGER.
*           On entry, NE specifies how many stations will be treated 
*           as eligible for consideration. 
*           Unchanged on exit.
*
*  NS       (input) INTEGER.
*           On entry, NS specifies how many stations have to be added 
*           to the current network. 1 <= NS <= NE.
*           Unchanged on exit.
*
*  AS       (work array) DOUBLE PRECISION, DIMENSION (NA,NA).
*
*  BS       (work array) DOUBLE PRECISION, DIMENSION (NA,NA).
*
*  CS       (work array) DOUBLE PRECISION, DIMENSION (NA,NA).
*
*  INV      (work array) DOUBLE PRECISION, DIMENSION (NA,NA).
*
*  IND      (work vector) INTEGER, DIMENSION (NA).
*
*  IERR     (output) INTEGER
*           the error indicator from DPOFA.  this will be
*           zero unless the matrix being factorized is
*           not positive definite.
*
*  W        (work vector) DOUBLE PRECISION, DIMENSION (NE).
*
*  WORK     (work vector) DOUBLE PRECISION, DIMENSION (LWORK).
*
*  LWORK    (integer) DIMENSION of the array WORK.
*
*  IWORK    (work vector) INTEGER, DIMENSION (8*NA).
*
*
*  =====================================================================
*
*     ..
*     .. Parameters ..
*     ..
*     .. Local Scalars ..
      INTEGER            I, NF, LF, LE, ILOEV, IUPEV, NEV, MI, NI, LFE
      DOUBLE PRECISION   SDET, ABSTOL, ALPHA, BETA, PROD
      CHARACTER*1        TRANSA,TRANSB,JOBZ,RANGE,UPLO
*     ..
*     .. Local Arrays ..
      INTEGER            IFAIL(1)
      DOUBLE PRECISION   DET(2), Z(1)
*     ..
*     .. External Functions ..
*     ..
*     .. External Subroutines ..
      EXTERNAL           PSUBM, DSCHUR, DSYEVX
*     ..
*     .. Intrinsic Functions ..
*     ..
*     .. Executable Statements ..
*
*
*     actual index sets of E and F

      NF=NA-NE
      DO 10 I=1,NF
         IND(I)=I
 10   CONTINUE

      LF=NF
      LFE=NF
      DO 20 I=1,NE
         IF(F(I).NE.0) THEN
            LFE=LFE+1
            LF=LF+1
            IND(LFE)=I+NF
         END IF
 20   CONTINUE
      LE=0
      DO 30 I=1,NE
         IF(E(I).NE.0) THEN
            LFE=LFE+1
            LE=LE+1
            IND(LFE)=I+NF
         END IF
 30   CONTINUE

      CALL PSUBM(A,LDA,NA,AS,LDAS,IND,LFE)

*     determining the Schur complement 
*     B[F,E]=A[E,E]-A[E,F]*A[F,F]^(-1)*A[F,E]
*     and the determinant of A[F,F]

      CALL DSCHUR(LFE,LF,AS,LDAS,BS,LDBS,SDET,IERR)

*     s-f greatest eigenvalues of B[F,E]

      JOBZ='N'
      RANGE='I'
      UPLO='U'
      ILOEV=LE-NS-NF+LF+1
      IUPEV=LE
      NEV=NS+NF-LF
      ABSTOL=1.0D-16
      CALL DSYEVX(JOBZ,RANGE,UPLO,LE,BS,LDBS,0.0D0,0.0D0,ILOEV,IUPEV,
     .            ABSTOL,NEV,W,Z,1,WORK,LWORK,IWORK,IFAIL,IERR)

*     product of the eigenvalues

      PROD=1.0
      DO 40 I=1,NEV
         PROD=PROD*W(I)
 40   CONTINUE     

*     new upper bound
      
      
      UPBND=SDET * PROD
      RETURN
      END
