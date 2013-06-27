module weights

contains 
   subroutine getWeights(dists,maxorder)
      double precision, intent(in) dists(:,:)
      double precision, intent(out) weights(size(A,1),size(A,2))
      integer i,j,m,n
      m=size(dists,1)
      n=size(dists,2)

      do i=2,m
         
!Advance the fabric evolution at a particular site
!Uses simple forward euler.. Maybe ok?
   subroutine solveJefferys(dt,epsdot,spin,ICs)
      double precision, intent(in) :: epsdot(3,3),spin(3,3)
      double precision, intent(inout) :: ICs(2,:)
      double precision :: n(3)

      n=
      
