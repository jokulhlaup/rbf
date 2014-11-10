module weights
   contains 
!   subroutine getWeights(dists,maxorder)
!      implicit none
!      double precision, intent(in) dists(:,:)
!      double precision, intent(out) weights(size(A,1),size(A,2))
!      integer i,j,m,n
!      m=size(dists,1)
!      n=size(dists,2)
!      end subroutine



!Advance the fabric evolution at a particular site
!Uses explicit RK4
!Only pass the first two rows of epsdot at spin.
!This consumes theta and phi
!also this assumes that the strain rate varies slow enough
!that we can take it as constant over one runge-kutta step.
!ie, this is first order accurate in spin,epsdot,
!but fourth order for (theta,phi)
   subroutine solveJefferys(dt,epsdot,spin,theta,phi,p)
   !!!Args
   !dt <- timestep
   !p <- number of xtals at single site. Needed for f2py.
   !epsdot <- strain rate tensor
   !spin <- spin tensor
   !theta,phi <- spherical coors [0:2pi], [0,pi] resp
      implicit none
      integer :: i,p
      double precision, intent(in) :: epsdot(3,3),spin(3,3)
      double precision, intent(inout) :: theta(p),phi(p)
      double precision :: n(3),dt,k1(3),k2(3),k3(3),k4(3)
      !loop over each xtal at site
      do i=1,p
         n(1)=sin(theta(i))*cos(phi(i))
         n(2)=sin(theta(i))*sin(phi(i))
         n(3)=cos(theta(i))
         k1=dt*rhs(spin,epsdot,n)
         k2=dt*rhs(spin,epsdot,n+0.5*k1)
         k3=dt*rhs(spin,epsdot,n+0.5*k2)
         k4=dt*rhs(spin,epsdot,n+0.5*k3)
         n=n+(k1+2*(k2+k3)+k4)/6
         !n=(matmul(spin,n)-( matmul(epsdot,n) - dot_product(matmul(epsdot,n),n)*n))*dt+n


         !now move back to spherical coordinates
         theta(i)=acos(n(3))
         phi(i)=asin(n(2)/sin(theta(i)))
         end do
      end subroutine 
  function rhs(spin,epsdot,n) result(r)
     implicit none
     double precision,intent(in) :: spin(3,3),epsdot(3,3),n(3)
     double precision r(3)
     r=(matmul(spin,n)-( matmul(epsdot,n) - dot_product(matmul(epsdot,n),n)*n))
     end function
end module 
