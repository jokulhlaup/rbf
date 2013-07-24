module jefferys
using ODE
export FabricPt,GlobalPars,solveJefferys,rk4,nRK4,rotC

 


#Modification of ODE4 from package ODE
function nRK4(f,ntimes,h,m,p)
  for i=1:ntimes
     p[:,i]=jefferys.rk4(f,h,m,p[:,i])
     end
  return p
  end

function rk4(f::Function,h::Float64,n::Int64,x::Array{Float64,1})
   for i=1:n
      k1=f(x)
      k2=f(x+k1*h/2)
      k3=f(x+k2*h/2)
      k4=f(x+k3*h/2)
      x+=(1/6)*h*(k1+2*k2+2*k3+k4)
      end
   return x
   end



##########################
##########Get viscosity###
##########################
type FabricPt
  p::Array{Float64,2} #[2,:] (theta,phi) angles
  coors::Array{Float64,1} #coors in space
  n::Int64 #number of xtals at site
  end

immutable GlobalPars
  dt::Float64 #timestep between velocity timesteps
  nrk::Int64 #Number of timesteps to be taken per dt for Jeffery's eqn by RK4
  hrk::Float64 #better be dt/nrk
  f::Function #Jeffery's eqn to supply to rk3
  end

#Replace this so its rotC(R)
function rotC(R)
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
  C = zeros(6,6)
  C[5,5]=1 
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
  p=nRK4(pars.f,fab.n,pars.hrk,pars.nrk,fab.p)
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
