module jefferys
using ODE
export Fabric,GlobalPars,solveJefferys,rk4,nRK4,rotC
##########################
##########Get viscosity###
##########################
abstract AbstractFabric<:Any

type Fabric{T<:Number,I<:Int}<:AbstractFabric
  coors::Array{T,2} #coors in space
  p::Array{T,3} #[2,:] (theta,phi) angles
  ngr::I #number of grains at site
  h::T
  ns::I #number of sites
  C::Array{T,3} #viscosity matrix
  #Fabric(coors,p,ngr,h,ns,C=zeros)=new(coors,p,ngr,h,ns,C)
  #stencil::Array{T,1}
 #Fabric(coors,p,ngr,ns,C)=new(coors,p,n,C,stencil)
  function Fabric(coors,p,ngr,ns,h,C=zeros(ns,6,6))
    size(coors)==(ns,3)?nothing:error("Dimension mismatch in coors")
    size(p)==(ns,ngr,3)?nothing:error("Dimension mismatch in p")
    size(C)==(ns,6,6) ?nothing:error("Dimension mismatch in C")
    return new(coors,p,ngr,h,ns,C)
    end
  end
immutable GlobalPars{T<:Number,I<:Int}
  dt::T #timestep between velocity timesteps
  nrk::I #Number of timesteps to be taken per dt for Jeffery's eqn by RK4
  hrk::T #better be dt/nrk
  f::Function #Jeffery's eqn to supply to rk3
  GlobalPars(dt,nrk,hrk,f)=new(dt,nrk,dt/nrk,f)
  end

#Modification of ODE4 from package ODE
function nRK4(f,ntimes,h,m,p)
  for i=1:ntimes
     p[i,:]=jefferys.rk4(f,h,m,p[i,:])
     end
  return p
  end

function rk4(f::Function,h::Float64,n::Int64,x::Array{Float64,1},vort,epsdot,theta,dt,m)
   for i=1:n
      k1=f(x)
      k2=f(x+k1*h/2)
      k3=f(x+k2*h/2)
      k4=f(x+k3*h/2)
      x+=(1/6)*h*(k1+2*k2+2*k3+k4)
      end
   return x
   end

function jefferysRHS(c,vort,epsdot,theta,dt,m)
  return (vort*c-epsdot*c+(c'*epsdot*c)*c)
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

#gets the rotation matrices
function getRotM(fab::T<:AbstractFabric)
  R=Array(Float64,fab.ngr,3,3)
  for i=1:fab.ngr
    A=[0 0 fab.p[i,1]
       0 0 -fab.p[i,2] 
       -fab.p[i,1] fab.p[i,2] 0]
    A2=[-fab.p[i,2]^2 fab.p[i,2]*fab[i,1] 0
        fab.p[i,2]*fab[i,1] fab.p[i,2]^2 0 
        0 0 fab.p[i,1]^2+fab.p[i,2]^2]
    R[i,:,:]=sin(acos(fab.p[i,3]))*A+(1-fab.p[i,3])*A2
    R[i,1,1]+=1;R[i,2,2]+=1;R[i,3,3]+=1
    end
  return R
  end
#this is the closure that returns fabEvolve!, which evolves the fabric based on the timestep.
function fabricHelper(pars::GlobalPars,f::Function,coors::Array{Float64,2},p::Array{Float64,2},ngrain::Int64,C::Array{Float64,3})
  #this is the function that actually does the rotation.  
  function fabEvolve!(fab::T<:AbstractFabric)
    fab.p=nRK4(fab.p,fab.ngr*fab.ns,fab.h,fab.m,fab.p)
    end
  return fabEvolve!
  end
(f,ntimes,h,m,p)


#Main driver routine to get the viscosity.
#refactor this by 
#at step zero, generate a closure equiv to fabEvolve! +params
#then at each timestep, mutate the closure. (can you do that
#without pushing a new copy onto the stack?)
function getVisc!(fab::T<:AbstractFabric,pars::GlobalPars,fabEvolve::Function)
  #advance the viscosity
  #get new theta
  fabEvolve(fab,pars)
  #get the rotation matrices for each theta pair
  R=getRotM(fab)
  #now rotate each grain's C matrix (C_ij= delta5,5)
  #Then invert it to get visc. yay.
  for i=1:fab.ngr
    C[i,:,:,]=rotC(R[i,:,:])
    end
  aC=mean(C,1) #check if right
  fab.visc=inv(aC)
  end



end #module
