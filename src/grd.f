      SUBROUTINE GRD(A,LDA,NA,NF,NE,NS,S,OPT,IND,AS,LDAS,BS,LDBS,INV,
     .               LDINV,V,W,IERR)
      IMPLICIT NONE
*
*  -- user defined driver routine --   
*     
*
*     .. Scalar Arguments ..
      INTEGER            NA, NF, NE, NS, IERR, LDA, LDAS, LDBS, LDINV
      DOUBLE PRECISION   OPT
*     ..
*     .. Array Arguments ..
      INTEGER            S(*), IND(*)
      DOUBLE PRECISION   A(LDA,*), AS(LDAS,*), BS(LDBS,*), INV(LDINV,*)
      DOUBLE PRECISION   V(*), W(*)
*     ..
*
*  Purpose
*  =======
*
*  GRD computes a Greedy solution for maximum entropy sampling
*  using the greedy algorithm described in Ko, Lee, Queyranne (1994),         
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
*           every feasible solution. 1 <= NF <= NA. 
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
*           to the current network. 1 <= NS <= NE.
*           Unchanged on exit.
*
*  S        (output) INTEGER array of DIMENSION (NE).
*           Contains The initial solution computed by the 
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
*
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
      INTEGER            I, J, K, CARDS, JMIN, IMAX, NI, INC, I0
      DOUBLE PRECISION   DET1, DETA, DETMAX, RCOND, ALPHA, BETA, PROD
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
      EXTERNAL           SUBDIN, PSUBV, DGEMV
*     ..
*     .. Intrinsic Functions ..
*     ..
*     .. Executable Statements ..
*
*
* Initializing 

      DETMAX=0.0D0
      DO 10 I=1,NE
         S(I)=0
 10   CONTINUE
      DO 11 I=1,NF
         IND(I)=I
 11   CONTINUE
      DO 12 I=NF+1,NA
         IND(I)=0
 12   CONTINUE
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

* Is the actual solution S feasible?

      IF (CARDS.LT.(NF+NS)) THEN 
         GOTO 22
      ELSE
         GOTO 30
      ENDIF
 22   CONTINUE

* Determinant and inverse INV of the actual infeasible solution
      
      DETMAX=0.0D0
      
      NI=K-1
      CALL SUBDIN(A,LDA,NA,AS,LDAS,IND,NI,BS,LDBS,INV,LDINV,IERR,DET)
      DETA=DET(1)*10**DET(2)
      
* Selecting the station I0 as candidate for augmenting the network
* Computation of the determinant of the covariance matrix
* of the resulting network using cholesky factorization  
      
      DO 23 I=1,NE
         IF (S(I) .EQ. 0) THEN 
            I0=NF+I
            CALL PSUBV(A,LDA,NA,V,IND,NI,I0)
            ALPHA=1.0D0
            BETA=0.0D0
            INC=1
            TRANS='N'
            CALL DGEMV(TRANS,NI,NI,ALPHA,INV,LDINV,V,INC,BETA,W,INC)
            PROD=DDOT(NI,V,INC,W,INC)
            DET1=DETA*(A(NF+I,NF+I)-PROD)

* Looking for the maximal determinant
* Storing the index for which the determinant attained its maximum

            
            IF (DET1 .GE. DETMAX) THEN
               DETMAX=DET1
               IMAX=I
            ENDIF
         ENDIF
 23   CONTINUE

* Augmenting the network selecting the station IMAX

      S(IMAX)=NF+IMAX
      GOTO 20


 30   CONTINUE
      OPT=DETMAX
      DO 31 I=1,NS
 31   CONTINUE
      RETURN
      END 


      
    

 







