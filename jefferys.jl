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


function halton(n,dim,base=nothing,skip=1e3)
  base=base==nothing?:base:[1,3,5,7,11,13,17][1:dim]
  f=1/base
  i=index
  for ndim=1:dim
    while idim>0
      xs=xs+f*(i%base[idim])
      i=floor(i/base[idim])
      f=f/base[idim]
    end
  return xs
##########################
##########Get viscosity###
##########################
type FabricPt
  coors::Array{Float64,1} #coors in space
  p::Union(Nothing,Array{Float64,2}) #[2,:] (theta,phi) angles
  n::Union(Int64,Nothing) #number of xtals at site
  C::Union(Nothing,Array{Float64,2}) #viscosity matrix
  stencil::Union(Nothing,Array{Int32,1})
  FabricPt(coors,p=nothing,n=nothing,C=nothing,stencil=nothing)=new(coors,p,n,C,stencil)
  end

immutable GlobalPars
  dt::Float64 #timestep between velocity timesteps
  nrk::Int64 #Number of timesteps to be taken per dt for Jeffery's eqn by RK4
  hrk::Float64 #better be dt/nrk
  f::Function #Jeffery's eqn to supply to rk3
  GlobalPars(dt,nrk,hrk,f)=new(dt,nrk,dt/nrk,f)
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
function getRotM(fab::FabricPt)
  R=Array(Float64,3,3,fab.n)
  for i=1:fab.n
    A=[0 0 fab.p[1,i]
       0 0 -fab.p[2,i] 
       -fab.p[1,i] fab.p[2,i] 0]
    A2=[-fab[2,i]^2 fab[2,i]*fab[1,i] 0
        fab[2,i]*fab[1,i] fab[2,i]^2 0 
        0 0 fab[1,i]^2+fab[2,i]^2]
    R[:,:,i]=sin(acos(fab.p[3,i]))*A+(1-fab.p[3,i])*A2
    R[1,1,i]+=1;R[2,2,i]+=1;R[3,3,i]+=1
    end
  return R
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
