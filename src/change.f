      SUBROUTINE CHANGE(A,LDA,NA,NF,NE,NS,S,OPT,IND,AS,LDAS,BS,LDBS,INV,
     .                  LDINV,V,W,IERR)
      IMPLICIT NONE
*
*  -- user defined driver routine --   
*     
*
*     .. Scalar Arguments ..
      INTEGER            NA, NF, NE, NS, IERR, LDA, LDAS, LDINV, LDBS
      DOUBLE PRECISION   OPT
*     ..
*     .. Array Arguments ..
      INTEGER            S(*), IND(*)
      DOUBLE PRECISION   A(LDA,*), AS(LDAS,*), BS(LDBS,*), INV(LDINV,*),
     .                   V(*), W(*)
*     ..
*
*  Purpose
*  =======
*
*  CHANGE performes interchange steps to improve a Greedy or Dual Greedy 
*  solution for maximum entropy sampling using the greedy algorithm 
*  described in Ko, Lee, Queyranne (1994),         
*  An exact algorithm for maximum entropy sampling, Operations
*  Research 43 (1995), 684-691.       
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
*  NF       (input) INTEGER.
*           On entry, NF specifies how many indices are forced into 
*           every feasible solution. 0 <= NF <= NA and NF+NS>=2. 
*           The stations related with NF have to be the first one in A.
*           Unchanged on exit.
*
*  NE       (input) INTEGER.
*           On entry, NE specifies how many stations will be treated 
*           as eligible for consideration. NE = NA - NF. The matrix A
*           stations related with NE have to be the last one in A.
*           Unchanged on exit.
*
*  NS       (input) INTEGER.
*           On entry, NS specifies how many stations have to be added 
*           to the current network. 1 <= NS <= NE and NF+NS>=2.
*           Unchanged on exit.
*
*  S        (output) INTEGER array of DIMENSION (NE).
*           Contains the initial solution computed by the 
*           greedy algorithm.
*
*  OPT      (output) DOUBLE PRECISION value.
*           The determinant of the initial solution computed by the 
*           greedy algorithm.
*
*  IND      (work array) INTEGER, DIMENSION (NA).
*
*  AS       (work array) DOUBLE PRECISION, DIMENSION (NA,NA).
*
*  BS       (work array) DOUBLE PRECISION, DIMENSION (NA,NA).
*
*  INV      (work array) DOUBLE PRECISION, DIMENSION (NA,NA).
*
*  V        (work array) DOUBLE PRECISION, DIMENSION (NA).

*  W        (work array) DOUBLE PRECISION, DIMENSION (NA).
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
      INTEGER            I, J, K, CARDS, IMAX, NI, INC, J0
      DOUBLE PRECISION   DET1, DETA, DETMAX, RCOND, ALPHA, BETA, PROD,
     .                   DETTEST
      CHARACTER*1        TRANS
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   DET(2)
*     ..
*     .. External Functions ..
      DOUBLE PRECISION DDOT
      EXTERNAL DDOT
*     ..
*     .. External Subroutines ..
      EXTERNAL           SUBDET, PIND, SUBDIN, PSUBV, DGEMV
*     ..
*     .. Intrinsic Functions ..
*     ..
*     .. Executable Statements ..
*
*
* Initializing 
 
      DO 11 I=1,NF
         IND(I)=I
 11   CONTINUE
      IMAX=1
 
* Computation of the actual index set IND and its cardinality

 20   CONTINUE
      K=NF+1
      CARDS=NF
      DO 21 I=1,NE
         IF (S(I).NE.0) THEN
            CARDS=CARDS+1
            IND(K)=S(I)
            K=K+1
         ENDIF
 21   CONTINUE
      NI=K-1

 30   CONTINUE

* Computation of the determinant of the covariance matrix 
* of the actual solution (DETMAX)
* Storing the index IMAX

      DO 31 I=1,NE
         IMAX=I
         CALL SUBDET(A,LDA,NA,AS,LDAS,IND,NI,BS,LDBS,IERR,DET)
         DETMAX=DET(1)*10**DET(2)

* Delete station IMAX from the current network provided IMAX 
* belongs to the network
* Computation of the new index set IND
* Computation of the determinant of the covariance matrix
* of the actual network using cholesky factorization 

         IF (S(I) .NE. 0) THEN
            S(I)=0
            CALL PIND(S,IND,NA,NE,NF,NI)
            
            CALL SUBDIN(A,LDA,NA,AS,LDAS,IND,NI,BS,LDBS,INV,LDINV,
     .           IERR,DET)
            DETA=DET(1)*10**DET(2)

* It is possible to improve the actual optimum DETMAX, 
* if we replace IMAX with J? 
* Computation of the determinant of the covariance matrix
* of the resulting network using cholesky factorization 
            
            DO 32 J=1,NE
               IF (S(J) .EQ. 0 .AND. J .NE. I) THEN
                  S(J)=NF+J
                  J0=NF+J
                  CALL PSUBV(A,LDA,NA,V,IND,NI,J0)
                  ALPHA=1.0D0
                  BETA=0.0D0
                  INC=1
                  TRANS='N'
                  CALL DGEMV(TRANS,NI,NI,ALPHA,INV,LDINV,V,INC,BETA,
     .                 W,INC)
                  PROD=DDOT(NI,V,INC,W,INC)
                  DET1=DETA*(A(NF+J,NF+J)-PROD)
                  
* Looking for the maximal determinant
* Storing the index for which the determinant attained its maximum
* Calculating the test determinant DETTEST avoids a infinite loop in case 
* of rounding errors occuring for (theoretically) identical determinants
                  
                  IF (DET1 .GT. DETMAX) THEN
                     CALL PIND(S,IND,NA,NE,NF,NI)
                     CALL SUBDET(A,LDA,NA,AS,LDAS,IND,NI,BS,LDBS,IERR,
     .                    DET)
                     DETTEST=DET(1)*10**DET(2)
                     IF (DETTEST .GT. DETMAX) THEN
                        DETMAX=DETTEST
                        IMAX=J
                     ENDIF
                  ENDIF
                  S(J)=0
               ENDIF
 32         CONTINUE

* Station IMAX belongs to the (new) network 
* Computation of the actual index set IND

            S(IMAX)=NF+IMAX
            CALL PIND(S,IND,NA,NE,NF,NI)
            
* If IMAX was changed leave the loop and start from beginning 
            
            IF (IMAX .NE. I) GOTO 30
         ENDIF
 31   CONTINUE
      
      OPT=DETMAX

      DO 41 I=1,NS
 41   CONTINUE
      END 


      
    

 







