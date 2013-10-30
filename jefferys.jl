module jefferys
using ODE
export Fabric,Fabric2,GlobalPars,AbstractFabric,solveJefferys,rk4,nRK4,rotC,jefferysRHS,fabricHelper
##########################
##########Get viscosity###
##########################
abstract AbstractFabric{T,I}<:Any
#most basic
type Fabric{T<:Number,I<:Int}<:AbstractFabric
  coors::Array{T,2} #coors in space
  p::Array{T,3} #[2,:] (theta,phi) angles
  ngr::I #number of grains at site
  h::T
  ns::I #number of sites
  C::Array{T,3} #viscosity matrix
  vort::Array{T,3} #vorticity
  epsdot::Array{T,3} #strain rate
  
  #Fabric(coors,p,ngr,h,ns,C=zeros)=new(coors,p,ngr,h,ns,C)
  #stencil::Array{T,1}
 #Fabric(coors,p,ngr,ns,C)=new(coors,p,n,C,stencil)
  function Fabric(coors,p,ngr,ns,h,C,vort,epsdot)
    size(coors)==(3,ns)?nothing:error("Dimension mismatch in 'coors'")
    size(p)==(3,ns,ngr)?nothing:error("Dimension mismatch in 'p'")
    size(C)==(6,6,ns)?nothing:error("Dimension mismatch in 'C'")
    size(vort)==(3,3,ns)?nothing:error("Dimension mismatch in 'vort'")
    size(epsdot)==(3,3,ns)?nothing:error("Dimension mismatch in 'epsdot'")
    return new(coors,p,ngr,h,ns,C,vort,epsdot)
    end
  end

type Fabric2{T<:Number,I<:Int}<:AbstractFabric
  coors::Array{T,2} #coors in space
  p::Array{T,3} #[2,:] (theta,phi) angles
  ngr::I #number of grains at site
  ns::I #number of sites
  h::T
  prob::T
  C::Array{T,3} #viscosity matrix
  vort::Array{T,3} #vorticity
  epsdot::Array{T,3} #strain rate
  #Fabric(coors,p,ngr,h,ns,C=zeros)=new(coors,p,ngr,h,ns,C)
  #stencil::Array{T,1}
 #Fabric(coors,p,ngr,ns,C)=new(coors,p,n,C,stencil)
  function Fabric2(coors,p,ngr,ns,h,prob,C,vort,epsdot)
    size(coors)==(3,ns)?nothing:error("Dimension mismatch in 'coors'")
    size(p)==(3,ns,ngr)?nothing:error("Dimension mismatch in 'p'")
    size(C)==(6,6,ns)?nothing:error("Dimension mismatch in 'C'")
    size(vort)==(3,3,ns)?nothing:error("Dimension mismatch in 'vort'")
    size(epsdot)==(3,3,ns)?nothing:error("Dimension mismatch in 'epsdot'")
    return new(coors,p,ngr,h,ns,prob,C,vort,epsdot)
    end
  end


immutable GlobalPars{T<:Number,I<:Int}
  dt::T #timestep between velocity timesteps
  nrk::I #Number of timesteps to be taken per dt for Jeffery's eqn by RK4
  hrk::T #better be dt/nrk
  f::Function #Jeffery's eqn to supply to rk3
  function GlobalPars(dt,nrk,f)
    return new(dt,nrk,dt/nrk,f)
    end
  end

#Modification of ODE4 from package ODE
function nRK4(f,ntimes::Int,h::Number,m::Number,dt,p,vort,epsdot)
  for i=1:ntimes
     p[:,i]=jefferys.rk4(f,h,m,p[:,i],vort[:,:],epsdot[:,:],dt,m)
     end
  return p
  end

function rk4(f::Function,h::Float64,n::Int64,x,vort,epsdot,dt,m)
   for i=1:n
      k1=f(x,vort,epsdot,dt)
      k2=f(x+k1*h/2,vort,epsdot,dt)
      k3=f(x+k2*h/2,vort,epsdot,dt)
      k4=f(x+k3*h/2,vort,epsdot,dt)
      x+=(1/6)*h*(k1+2*k2+2*k3+k4)
      end
   return x
   end

function jefferysRHS(c,vort,epsdot,dt)
  return (epsdot*c-(c'*epsdot*c)[1]*c)*dt
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
function getRotM(fab::AbstractFabric)
  R=Array(Float64,fab.ngr,3,3)
  for i=1:fab.ngr*fab.ns
    A=[0 0 fab.p[1,i]
       0 0 -fab.p[2,i] 
       -fab.p[1,i] fab.p[2,i] 0]
    A2=[-fab.p[2,i]^2 fab.p[2,i]*fab[1,i] 0
        fab.p[2,i]*fab[1,i] fab.p[2,i]^2 0 
        0 0 fab.p[1,i]^2+fab.p[2,i]^2]
    R[:,:,i]=sin(acos(fab.p[3,i]))*A+(1-fab.p[3,i])*A2
    R[1,1,i]+=1;R[2,2,i]+=1;R[3,3,i]+=1
    end
  return R
  end
#this is the closure that returns fabEvolve!, which evolves the fabric based on the timestep.
function fabricHelper(pars::GlobalPars,fab::AbstractFabric,f::Function)
  #this is the function that actually does the rotation.  
  function fabEvolve!(pars::GlobalPars,fab::Fabric,f) #jefferys equation
    for i=1:fab.ns
      fab.p[:,i,:]=nRK4(f,fab.ngr,fab.h,pars.nrk,pars.dt,fab.p[:,i,:],fab.vort[:,:,i],fab.epsdot[:,:,i])
      end
      return fab.p
    end

  function fabEvolve!(pars::GlobalPars,fab::Fabric2,f) #jefferys equation
    for i=1:fab.ns
      fab.p[:,i,:]=nRK4(f,fab.ngr,fab.h,pars.nrk,pars.dt,fab.p[:,i,:],fab.vort[:,:,i],fab.epsdot[:,:,i])
      end
      dynRextal!(fab::Fabric2)
      return fab.p
  end
  function dynRextal!(fab::Fabric2)
    change=rand(fab.ns*fab.ngr)
    for i=1:fab.ns*fab.ngr
      p=rand()
      if p<fab.prob
        x=rand(3)-0.5
        fab.p[3*i-2:3*i]=x/norm(x)
        end
      end
    end
    return fabEvolve!
  end


#Main driver routine to get the viscosity.
#refactor this by 
#at step zero, generate a closure equiv to fabEvolve! +params
#then at each timestep, mutate the closure. (can you do that
#without pushing a new copy onto the stack?)
function getVisc!(fab::AbstractFabric,pars::GlobalPars,fabEvolve::Function)
  #advance the viscosity
  #get new theta
  fabEvolve!(fab,pars)
  #get the rotation matrices for each theta pair
  R=getRotM(fab)
  #now rotate each grain's C matrix (C_ij= delta5,5)
  #Then invert it to get visc. yay.
  for i=1:fab.ngr
    C[:,:,i]=rotC(R[:,:,i])
    end
  aC=mean(C,1) #check if right
  fab.visc=inv(aC)
  end

end #module
