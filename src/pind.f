      SUBROUTINE PIND(S,IND,NA,NE,NF,NI)
      IMPLICIT NONE
*
*  -- user defined driver routine --   
*     
*
*     .. Scalar Arguments ..
      INTEGER            NA, NF, NE, NI
*     ..
*     .. Array Arguments ..
      INTEGER            S(NE), IND(NA)
*     ..
*
*  Purpose
*  =======
*
*  PIND computes proceeding from the actual solution only the indices
*  belonging to this solution.
*
*  Arguments
*  =========
*
*
*  S        (input) INTEGER vector of DIMENSION (NE).
*           On entry, S contains the initial solution computed by 
*           the greedy algorithm.
*           Unchanged on exit.
*
*  IND      (output) INTEGER vector, DIMENSION (NA).
*           On exit, the first NI arguments of the IND store
*           the indices contained in the actual solution.  
*
*  NA       (input) INTEGER.
*           On entry, NA specifies the maximal dimension of IND.
*           Unchanged on exit.
*
*  NE       (input) INTEGER.
*           On entry, NE specifies the dimension of S. NE = NA-NF. 
*           Unchanged on exit.
*
*  NF       (input) INTEGER.
*           On entry, NF specifies how many indices are forced into 
*           every feasible solution. 0 <= NS <= NA. 
*           Unchanged on exit.
*
*  NI       (output) INTEGER.
*           On exit, NI specifies the cardinality of the actual solution.
*
*
*  =====================================================================
*
*     ..
*     .. Parameters ..
*     ..
*     .. Local Scalars ..
      INTEGER            I, J, K
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
      K=NF+1
      DO 10 I=1,NE
         IF (S(I).NE.0) THEN
            IND(K)=S(I)
            K=K+1
         ENDIF
 10   CONTINUE
      NI=K-1
      END
