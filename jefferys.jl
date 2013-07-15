module jefferys
using ODE
export solveJefferys,rk4,nRK4
#Modification of ODE4 from package ODE

function nRK4(f,ntimes,h,m,theta)
  for i=1:ntimes
     #(tout,theta[:,i])=ode45(f,1,theta[:,i])
     theta[:,i]=jefferys.rk4(f,h,m,theta[:,i])
     end
     return theta
  end

function rk4(f::Function,h::Float64,n::Int64,x::Array{Float64,1})
   for i=1:n
      k1=f(x)
      k2=f(x+k1*h/2)
      k3=f(x+k2*h/2)
      k4=f(x+k3*h/2)
      x+=1./6.*h*(k1+2*k2+2*k3+k4)
      end
   return x
   end

function jefferysLHS(vort,epsdot,theta,dt,m)
   n=[sin(theta[1])*cos(theta[2]),sin(theta[1])*sin(theta[2]),cos(theta[1])] 
   return((vort*n)[2:3])-((epsdot*n)[2:3]-(n'*epsdot*n)[1]*n[2:3])
   end
#Using ODE
function solveJefferys{T}(vort::Array{Float64,2},epsdot::Array{Float64,2},theta::Array{Float64,1},dt::Float64,m::Int64)
  n=theta->[sin(theta[1])*cos(theta[2]),sin(theta[1])*sin(theta[2]),cos(theta[1])]
  f=theta->((vort*n(theta))[2:3])-((epsdot*n(theta))[2:3]-(n(theta)'*epsdot*n(theta))[1]*n(theta)[2:3])
#  f=theta0->jefferyLHS(vort,epsdot,theta0,dt,m)
  jefferys.rk4(f,dt/m,m,theta)
  return theta
  end

end


<<<<<<< HEAD
function rotC(C,R)
  # form the K matrix (based on Bowers 'Applied Mechanics of Solids', Chapter 3)
  K1 = [ R[1,1]^2 R[1,2]^2 R[1,3]^2 ; 
         R[2,1]^2 R[2,2]^2 R[2,3]^2 ; 
         R[3,1]^2 R[3,2]^2 R[3,3]^2 ] ;
  K2 = [ R[1,2]*R[1,3] R[1,3]*R[1,1] R[1,1]*R[1,2] ; 
         R[2,2]*R[2,3] R[2,3]*R[2,1] R[2,1]*R[2,2] ;
         R[3,2]*R[3,3] R[3,3]*R[3,1] R[3,1]*R[3,2] ] ;
  
  K3 = [ R[2,1]*R[3,1] R[2,2]*R[3,2] R[2,3]*R[3,3] ; 
         R[3,1]*R[1,1] R[3,2]*R[1,2] R[3,3]*R[1,3] ; 
         R[1,1]*R[2,1] R[1,2]*R[2,2] R[1,3]*R[2,3] ] ;
  
  K4 = [ R[2,2]*R[3,3]+R[2,3]*R[3,2] R[2,3]*R[3,1]+R[2,1]*R[3,3] R[2,1]*R[3,2]+R[2,2]*R[3,1] ; 
         R[3,2]*R[1,3]+R[3,3]*R[1,2] R[3,3]*R[1,1]+R[3,1]*R[1,3] R[3,1]*R[1,2]+R[3,2]*R[1,1] ;       
         R[1,2]*R[2,3]+R[1,3]*R[2,2] R[1,3]*R[2,1]+R[1,1]*R[2,3] R[1,1]*R[2,2]+R[1,2]*R[2,1]] ; 
  K = [ K1  2*K2 ; 
        K3   K4   ] ;
  
  CR = K * C * transpose(K) ;
  
  end
=======

##########################
##########Get viscosity###
######################### #
>>>>>>> 211d7f6b9dffde6a871a146b01d0ad1d7c623b92
