module jefferys
using ODE
export Fabric,Fabric2,FabricNGG,genrFT,GlobalPars,AbstractFabric,solveJefferys,rk4,nRK4,rotC,jefferysRHS,fabricHelper
##########################
##########Get viscosity###
##########################
abstract AbstractFabric{T,I}<:Any
#most basic

#type generating function
function  genrFT(name,body)
  eval(
  quote
    type $name{T<:Number,I<:Int}<:AbstractFabric
      coors::Array{T,2} #coors in space
      p::Array{T,3} #[2,:] (theta,phi) angles
      ngr::I #number of grains at site
      h::T
      ns::I #number of sites
      C::Array{T,3} #viscosity matrix
      vort::Array{T,3} #vorticity
      epsdot::Array{T,3} #strain rate
      
      $body
      end
      end)
  end
    #Fabric(coors,p,ngr,h,ns,C=zeros)=new(coors,p,ngr,h,ns,C)
    #stencil::Array{T,1}
   #Fabric(coors,p,ngr,ns,C)=new(coors,p,n,C,stencil)
genrFT(:(Fabric),:(begin
  function Fabric(coors,p,ngr,ns,h,C,vort,epsdot)
    size(coors)==(3,ns)?nothing:error("Dimension mismatch in 'coors'")
    size(p)==(3,ngr,ns)?nothing:error("Dimension mismatch in 'p'")
    size(C)==(6,6,ns)?nothing:error("Dimension mismatch in 'C'")
    size(vort)==(3,3,ns)?nothing:error("Dimension mismatch in 'vort'")
    size(epsdot)==(3,3,ns)?nothing:error("Dimension mismatch in 'epsdot'")
    return new(coors,p,ngr,h,ns,C,vort,epsdot)
    end
  end))

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
    size(p)==(3,ngr,ns)?nothing:error("Dimension mismatch in 'p'")
    size(C)==(6,6,ns)?nothing:error("Dimension mismatch in 'C'")
    size(vort)==(3,3,ns)?nothing:error("Dimension mismatch in 'vort'")
    size(epsdot)==(3,3,ns)?nothing:error("Dimension mismatch in 'epsdot'")
    return new(coors,p,ngr,h,ns,prob,C,vort,epsdot)
    end
  end

type FabricNGG2{T<:Number,I<:Int}<:AbstractFabric
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
  nbrs::Array{I,2} #Associates nbrs. Hopefully transitive.
  #nbrs[:,i] is the nbrs of the i'th grain, where i in [1:ngr*ns]
  r::Array{T,1} #radius of grains.
  #Probably don't want to mix grains from different sites

 #Fabric(coors,p,ngr,ns,C)=new(coors,p,n,C,stencil)
  function FabricNGG2(coors,p,ngr,ns,h,prob,C,vort,epsdot,nbrs,r)
    size(coors)==(3,ns)?nothing:error("Dimension mismatch in 'coors'")
    size(p)==(3,ngr,ns)?nothing:error("Dimension mismatch in 'p'")
    size(C)==(6,6,ns)?nothing:error("Dimension mismatch in 'C'")
    size(vort)==(3,3,ns)?nothing:error("Dimension mismatch in 'vort'")
    size(epsdot)==(3,3,ns)?nothing:error("Dimension mismatch in 'epsdot'")
    return new(coors,p,ngr,h,ns,C,vort,epsdot,nbrs,r)
    end
  end



genrFT(:(FabricNGG),:(begin
  nbrs::Array{I,2} #Associates nbrs. Hopefully transitive.
  #nbrs[:,i] is the nbrs of the i'th grain, where i in [1:ngr*ns]
  r::Array{T,1} #radius of grains.
  #Probably don't want to mix grains from different sites
  end))

function consFabricNGG(coors,p,ngr,ns,h,C,vort,epsdot,nn,av_radius)
  nbrs=makeRandomNbrs(ns,ngr,nn)
  r=2*rand(ngr*ns)*av_radius

  size(coors)==(3,ns)?nothing:error("Dimension mismatch in 'coors'")
  size(p)==(3,ngr,ns)?nothing:error("Dimension mismatch in 'p'")
  size(C)==(6,6,ns)?nothing:error("Dimension mismatch in 'C'")
  size(vort)==(3,3,ns)?nothing:error("Dimension mismatch in 'vort'")
  size(epsdot)==(3,3,ns)?nothing:error("Dimension mismatch in 'epsdot'")
  length(nbrs)==ngr*ns*nn?nothing:error("Dimension mismatch in 'nbrs'")
  return FabricNGG{Float64,Int64}(coors,p,ngr,h,ns,C,vort,epsdot,nbrs,r)
  end

######################
############
#####################
module Ngg
type Nbrs
  nbrs::Array{Int64,2}
  end

function advanceRadius(this,rs,grmob,dt)
  dr=nggRate(this,rs,grmob,dt)
  new=this+dr
  if new<0
    return 0
    else (return dr)
    end

  #get the velocity for ngg between two grains
  #Important: OUTWARD from r1, don't fuck that up.
  function nggVelocity(r1,r2,grmob)
    mc=(1./r2-1./r1)./2
    return mc*grmob
    end
  #returns radius dt
  function nggRate(this,rs,grmob,dt)
    pa=propAreas(rs)
    v=nggVelocity(this,rs,grmob)
    dV=2/3*pi*(this*this*this-dot(pa,(this-v*dt).^3))
    if dV<0
      return -((abs(dV))^3)
      else
      dV^(1/3)
      end
    end
  end 

function makeRandomNbrs(ns::Int,ngr::Int,nn::Int)
  nbrs=Nbrs(Array(Int,nn,ns*ngr))
  for i=

function deleteNbr(this::Int)
  

function setConjugateNbrs(this::Int,nbr::Array{Float64,1},nbrs::Nbrs,nn,repl::Bool=false)
  for i=1:nn
    if nbr[i]==0
      nbr[i]=this
      return nbr
      end
    end
  i=Utils.randir(1,nn)
  nbrs[i]=this
  end




function makeRandomNbrs(ns::Int,ngr::Int,nn::Int)
nbrs=Array(Int,nn,ns*ngr)
for i=1:ns
  for j=1:ngr
    for k=1:nn
      nbrs[k,(i-1)*ngr+j]=Utils.diffrandi((i-1)*ngr+j,(i-1)*ngr,i*ngr)
      end
    end
  end
  return nbrs
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
  K4 = [ R[2,2]*R[3,3]+R[2,3]*R[3,2] R[2,3]*R[3,1]+R[2,1]*R[3,3] 
         R[2,1]*R[3,2]+R[2,2]*R[3,1] ; 
         R[3,2]*R[1,3]+R[3,3]*R[1,2] R[3,3]*R[1,1]+R[3,1]*R[1,3] 
         R[3,1]*R[1,2]+R[3,2]*R[1,1] ;       
         R[1,2].*R[2,3]+R[1,3].*R[2,2] R[1,3].*R[2,1]+
         R[1,1].*R[2,3] R[1,1].*R[2,2]+R[1,2].*R[2,1]] ; 
  K = [ K1  2*K2 ; 
        K3   K4   ] ;
  C = zeros(6,6)
  C[5,5]=1 
  C = K * C * transpose(K) 
  end

#this is the closure that returns fabEvolve!,
#which evolves the fabric based on the timestep.
function fabricHelper(pars::GlobalPars,fab::AbstractFabric,f::Function)
  #this is the function that actually does the rotation.  
  function fabEvolve!(pars::GlobalPars,fab::Fabric,f) #jefferys equation
    for i=1:fab.ns*fab.ngr
      fab.p[:,i]=rk4(f,fab.ngr,fab.h,pars.nrk,pars.dt,fab.p[:,i,:],
          fab.vort[:,:,i],fab.epsdot[:,:,i])
      end
      return fab.p
    end
  
  #function fabEvolve!(pars::GlobalPars,fab::Fabric3,f)
  function fabEvolve!(pars::GlobalPars,fab::FabricNGG,f)
    for i=1:fab.ns*fab.ngr
      fab.p[:,i]=rk44(f,fab.ngr,fab.h,pars.nrk,pars.dt,fab.p[:,i],
        fab.vort[:,i],fab.epsdot[:,i])
      r[i]=advanceRadius(r[i],r[nbrs[:,i]],fab.grmob,fab.dt)
      end
    end
  function fabEvolve!(pars::GlobalPars,fab::Fabric2,f) #jefferys equation
    for i=1:fab.ns
      fab.p[:,:,i]=nRK4(f,fab.ngr,fab.h,pars.nrk,pars.dt,fab.p[:,:,i],
          fab.vort[:,:,i],fab.epsdot[:,:,i])
      end
      fab.p=dynRextal!(fab::Fabric2)
      return fab.p
    end
   return fabEvolve!
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

#finds the effective stress on grains at one site.
function localSigmaEff(p,sigma,ngr)
  G=zeros(6)
  #first find local geometric tensors
  for i=1:ngr
    g[6*i-5:6*i]=localGeomTensor(p[3*i-2:3*i],sigma)
    G+=g[6*i-5:6*i]
    end
  G=G/ngr
  for i=1:ngr
    g[6*i-5:6*i]=g[6*i-5:6*i]./G.*sigma #g= local sigma now
    sigmaE[i]=sqrt(1/3*secondInv(g[6*i-5:6*i])) #be sure to convert
    #back to effective stress from the second invariant.
    end
  return sigmaE 
  end
#probability that a crystal recrystallizes
function probDRx(fab::AbstractFabric)
  for i=1:ns
    sigmaE=localSigmaEff(fab.p[(i-1)*fab.ngr+1:i*fab.ngr],
        sigma[:,(i-1)*6+1:i*6],fab.ngr)
    A=expFactor()
    end
  end
#Find the local geometric tensor G (Azuma 1996)
#p::3 x ngr array
#g::6 x ngr array (symmetric)
function localGeomTensor(p,sigma)
  m=dot(sigma,c) #read: T
  m=cross(c,T) #read: n
  m=m/norm(m) #read: n
  m=cross(m,c)
  return m/norm(m)
  end
##############################  
#stuff for normal grain growth
#this gets the area proportions
function propAreas(rs)
  rs2=rs.*rs
  s=sum(rs2)
  if s==0
   return 1/length(rs)
   else 
     return rs2/sum(rs2)
   end
  end
#gets the rotation matrices
function getRotMHelper(Array{T,1},A::Array{Float64,2},
    A2=Array{Float64,2},R::Array{Float64,2})
  A=[0 0 p[1]
     0 0 p[2] 
     -p[1] p[2] 0]
  A2=[-p[2]^2 p[2]*fab[1] 0
      p[2]*p[1] p[2]^2 0 
      0 0 p[1]^2+p[2]^2]
  R[:,:]=sin(acos(p[3]))*A+(1-p[3])*A2
  R[1,1]+=1;R[2,2]+=1;R[3,3]+=1
  return R
  end

function getC(fab::AbstractFabric,A::Array{Float64,2},
    A2=Array{Float64,2},R::Array{Float64,2})
  for i=1:fab.ns*fab.ngr
    getRotMHelper(fab.p[3*i-2:3*i],A,A2,R)
    end
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
  R=Array(Float64,3,3)
  A=Array(Float64,3,3)
  A2=Array(Float64,3,3)
  #get the rotation matrices for each theta pair
  getC!(fab,A,A2,R)
  #now rotate each grain's C matrix (C_ij= delta5,5)
  #Then invert it to get visc. yay.
  aC=mean(C,3) #check if right
  fab.visc=inv(aC)
  end
end #module
