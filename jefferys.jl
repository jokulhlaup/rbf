module jefferys
using ODE
export solveJefferys,rk4,nRK4,rotC

 


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



##########################
##########Get viscosity###
##########################
type FabricPt
  theta::Array{Float64,2} #[2,:] (theta,phi) angles
  coors::Array{Float64,1} #coors in space
  n::Int64 #number of xtals at site
  end

type GlobalPars
  dt::Float64 #timestep between velocity timesteps
  nrk::Int64 #Number of timesteps to be taken per dt for Jeffery's eqn by RK4
  hrk::Float64 #better be dt/nrk
  f::Function #Jeffery's eqn to supply to rk3
  end

#Replace this so its rotC(R)
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
         R[1,2].*R[2,3]+R[1,3].*R[2,2] R[1,3].*R[2,1]+R[1,1].*R[2,3] R[1,1].*R[2,2]+R[1,2].*R[2,1]] ; 
  K = [ K1  2*K2 ; 
        K3   K4   ] ;
  C = K * C * transpose(K) 
end

function getRotM(fab::FabricPt)
  c1=cos(fab.theta[1])
  c2=cos(fab.theta[2])
  s1=sin(-fab.theta[1])
  s2=sin(-fab.theta[2])
  for i=1:fab.n 
    rotM[:,:,i]=
      [c2*c1 -s2 c2*s1;
       s2*c1 c2 s2*s1;
       -s1 0 c1]
    end
    return rotM
  end

function fabEvolve!(fab::FabricPt,pars::GlobalPars)
  #Move this to the outside as field of GlobalPars
  n=theta->[sin(theta[1])*cos(theta[2]),sin(theta[1])*sin(theta[2]),cos(theta[1])]
  f=theta->((vort*n(theta))[2:3])-((epsdot*n(theta))[2:3]-(n(theta)'*epsdot*n(theta))[1]*n(theta)[2:3])
  fab.theta=nRK4(f,fab.n,pars.hrk,pars.nrk,fab.theta)
  fab.theta[1,:]=fab.theta[1,:]%(2*pi)
  fab.theta[2,:]=fab.theta[2,:]%pi
  end 
#Main driver routine to get the viscosity.
function getVisc!(fab::FabricPt,pars::GlobalPars)
  #advance the viscosity
  #get new theta
  fabEvolve!(fab,pars)
  #get the rotation matrices for each theta pair
  R=getRotM(fab)
  #now rotate each crystal's C matrix (C_ij= delta5,5)
  #Then invert it to get visc. yay.
  for i=1:fab.n
    C[:,:,i]=rotC(R[:,:,i])
    end
    aC=mean(C,3)
    fab.visc=inv(aC)
  end



end #module
