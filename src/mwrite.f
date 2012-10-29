      SUBROUTINE mwrite(A,m,n)
      IMPLICIT NONE
      integer m,n
      double precision A(m,n)
      integer i,j

      do 10 i=1,m
         do 20 j=1,n
            write (*,*) i,",",j,":",A(i,j)
 20      continue
 10   continue
      end
