      SUBROUTINE GRDDL(A,LDA,NA,NF,NE,NS,S,OPT,IND,AS,LDAS,BS,LDBS,IERR)
      IMPLICIT NONE
*
*  -- user defined driver routine --   
*     
*
*     .. Scalar Arguments ..
      INTEGER            NA, NF, NE, NS, IERR, LDA, LDAS, LDBS
      DOUBLE PRECISION   OPT
*     ..
*     .. Array Arguments ..
      INTEGER            S(*), IND(*)
      DOUBLE PRECISION   A(LDA,*), AS(LDAS,*), BS(LDBS,*)
*     ..
*
*  Purpose
*  =======
*
*  GRDDL computes a Greedy solution for maximum entropy sampling
*  using the dual greedy algorithm described in 
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
*  NF       (input) INTEGER.
*           On entry, NF specifies how many indices are forced into 
*           every feasible solution. 0 <= NS <= NA. 
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
*           to the current network. 0 <= NS <= NE.
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
      INTEGER            I, J, K, CARDS, JMIN, JMAX, NI
      DOUBLE PRECISION   DET1, DETMAX
      CHARACTER*1        TRANS
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   DET(2)
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
* Initializing 

      DETMAX=0.0D0
      DO 10 I=1,NE
         S(I)=I+NF
 10   CONTINUE
      DO 11 I=1,NF
         IND(I)=I
 11   CONTINUE
      DO 12 I=NF+1,NA
         IND(I)=0
 12      CONTINUE
      JMIN=1

* Computation of the cardinality of the actual solution

 20   CONTINUE
      CARDS=NF
      DO 21 I=1,NE
         IF (S(I).NE.0) THEN
            CARDS=CARDS+1
         ENDIF
 21   CONTINUE

* Is the actual solution S feasible?

      IF (CARDS.GT.(NF+NS)) THEN 
         GOTO 22
      ELSE
         GOTO 30
      ENDIF
 22   CONTINUE

* Computation of the actual index set IND

      K=NF+1
      DO 221 I=JMIN,NE
         IF (S(I).NE.0 .AND. I.NE.JMIN) THEN 
            IND(K)=S(I)
            K=K+1
         ENDIF
 221  CONTINUE

* Computation of the determinant of the covariance matrix 
* of the actual infeasible solution

      NI=K-1
      CALL SUBDET(A,LDA,NA,AS,LDAS,IND,NI,BS,LDBS,IERR,DET)
      DETMAX=DET(1)*10.0**DET(2)
      JMAX=JMIN

* Selecting the station J as candidate to reduce the network
* Computation of the determinant of the covariance matrix
* of the resulting network using cholesky factorization 

      DO 23 J=JMIN+1,NE
         IF (S(J).NE.0) THEN
            K=NF+1
            DO 222 I=JMIN,NE
               IF (S(I).NE.0 .AND. I.NE.J) THEN 
                  IND(K)=S(I)
                  K=K+1
               ENDIF
 222        CONTINUE
            NI=K-1
            CALL SUBDET(A,LDA,NA,AS,LDAS,IND,NI,BS,LDBS,IERR,DET)
            DET1=DET(1)*10.0**DET(2)

* Looking for the maximal determinant
* Storing the index for which the determinant attained its maximum

            IF (DET1.GT.DETMAX) THEN 
               DETMAX=DET1
               JMAX=J
            ENDIF
         ENDIF
 23   CONTINUE

* Deleting the station JMAX from the current network

      S(JMAX)=0

* Determining the station with the smallest index in the current network   

      IF (JMIN.EQ.JMAX) THEN
         DO 24 J=JMAX+1,NE
            IF (S(J).NE.0) THEN
               GOTO 25
            ENDIF
 24      CONTINUE
 25      JMIN=J         
      ENDIF

      GOTO 20

 30   CONTINUE
      
      OPT=DETMAX
      DO 31 I=1,NS
 31   CONTINUE
      RETURN
      END 
