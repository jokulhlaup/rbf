module jefferys
using ODE,Utils, Debug
export consFabricNGG,Fabric,Fabric2,FabricNGG,genrFT,makeRandomNbrs!,fabEv!,advanceRadius,GlobalPars,AbstractFabric,solveJefferys,rk4,nRK4,rotC,jefferysRHS,fabricHelper
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
      sigma::Array{T,3} #strain rate
      
      $body
      end
      end)
  end
    #Fabric(coors,p,ngr,h,ns,C=zeros)=new(coors,p,ngr,h,ns,C)
    #stencil::Array{T,1}
   #Fabric(coors,p,ngr,ns,C)=new(coors,p,n,C,stencil)
genrFT(:(Fabric),:(begin
  function Fabric(coors,p,ngr,ns,h,C,vort,sigma)
    size(coors)==(3,ns)?nothing:error("Dimension mismatch in 'coors'")
    size(p)==(3,ngr,ns)?nothing:error("Dimension mismatch in 'p'")
    size(C)==(6,6,ns)?nothing:error("Dimension mismatch in 'C'")
    size(vort)==(3,3,ns)?nothing:error("Dimension mismatch in 'vort'")
    size(sigma)==(3,3,ns)?nothing:error("Dimension mismatch in 'sigma'")
    return new(coors,p,ngr,h,ns,C,vort,sigma)
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
  sigma::Array{T,3} #strain rate
  #Fabric(coors,p,ngr,h,ns,C=zeros)=new(coors,p,ngr,h,ns,C)
  #stencil::Array{T,1}
 #Fabric(coors,p,ngr,ns,C)=new(coors,p,n,C,stencil)
  function Fabric2(coors,p,ngr,ns,h,prob,C,vort,sigma)
    size(coors)==(3,ns)?nothing:error("Dimension mismatch in 'coors'")
    size(p)==(3,ngr,ns)?nothing:error("Dimension mismatch in 'p'")
    size(C)==(6,6,ns)?nothing:error("Dimension mismatch in 'C'")
    size(vort)==(3,3,ns)?nothing:error("Dimension mismatch in 'vort'")
    size(sigma)==(3,3,ns)?nothing:error("Dimension mismatch in 'sigma'")
    return new(coors,p,ngr,h,ns,prob,C,vort,sigma)
    end
  end

type FabricNGG2{T<:Number,I<:Int}<:AbstractFabric
  coors::Array{T,2} #coors in space
  p::Array{T,3} #[2,:] (theta,phi) angles
  ngr::I #number of grains at site
  ns::I #number of sites
  h::T
  C::Array{T,3} #viscosity matrix
  vort::Array{T,3} #vorticity
  sigma::Array{T,3} #strain rate
  #Fabric(coors,p,ngr,h,ns,C=zeros)=new(coors,p,ngr,h,ns,C)
  #stencil::Array{T,1}
  nn::I
  nbrs::Array{I,2} #Associates nbrs. Hopefully transitive.
  #nbrs[:,i] is the nbrs of the i'th grain, where i in [1:ngr*ns]
  r::Array{T,1} #radius of grains.
  #Probably don't want to mix grains from different sites

 #Fabric(coors,p,ngr,ns,C)=new(coors,p,n,C,stencil)
  function FabricNGG2(coors,p,ngr,ns,h,prob,C,vort,sigma,nbrs,r)
    size(coors)==(3,ns)?nothing:error("Dimension mismatch in 'coors'")
    size(p)==(3,ngr,ns)?nothing:error("Dimension mismatch in 'p'")
    size(C)==(6,6,ns)?nothing:error("Dimension mismatch in 'C'")
    size(vort)==(3,3,ns)?nothing:error("Dimension mismatch in 'vort'")
    size(sigma)==(3,3,ns)?nothing:error("Dimension mismatch in 'sigma'")
    return new(coors,p,ngr,h,ns,C,vort,sigma,nn,nbrs,r)
    end
  end



genrFT(:(FabricNGG),:(begin
  nn::I
  nbrs::Array{Bool,3}
  areas::Array{T,3} #Associates nbrs. Hopefully transitive.
  #nbrs[:,i] is the nbrs of the i'th grain, where i in [1:ngr*ns]
  r::Array{T,2} #radius of grains.
  grmob::T
  #Probably don't want to mix grains from different sites
  pr_nuc::T
  nuc_vol::T
  str::Array{T,2}
  end))

function consFabricNGG(coors,p,ngr,ns,h,C,vort,sigma,nn,av_radius)
  nbrs=makeSymNbrs(ns,ngr,0.1)
  r=2*rand(ngr,ns)*av_radius
  
  areas=zeros(ngr,ngr,ns)
  for i=1:ns
    areas[:,:,i]=initAreas(r[:,i],nbrs[:,:,i])
    end
  grmob=1.0
  pr_nuc=0.1
  nuc_vol=0.1
  str=zeros(ngr,ns)
  size(coors)==(3,ns)?nothing:error("Dimension mismatch in 'coors'")
  size(p)==(3,ngr,ns)?nothing:error("Dimension mismatch in 'p'")
  size(C)==(6,6,ns)?nothing:error("Dimension mismatch in 'C'")
  size(vort)==(3,3,ns)?nothing:error("Dimension mismatch in 'vort'")
  size(sigma)==(3,3,ns)?nothing:error("Dimension mismatch in 'sigma'")
  return FabricNGG{Float64,Int64}(coors,p,ngr,h,ns,C,vort,sigma,nn,nbrs,areas,r,grmob,pr_nuc,nuc_vol,str)
  end

######################
############
#####################
##########################################
#function advanceRadii(rs,nbrs,grmob,dt,sigma,ngr,areas,p,pr_nucleation,nuc_vol)
function advanceRadii(fab::FabricNGG,k,dt)
  
  #nucleates a new grain at position (i,j,k)
  function nucleateGrain!(fab,i,k,vol)
    jd=binBoolInd(fab.areas[:,i,k],>,ngr)
    vol[jd]-=fab.nuc_vol
    vol[i]=fab.nuc_vol
    if vol[jd]<0
      vol[i]+=vol[jd]
      vol[jd]=0
      end
    azimuth=rand()*2*pi  
    fab.p[:,i,k]=[sin(pi/4)*cos(azimuth),sin(pi/4)*sin(azimuth),cos(pi/4)]#rand(3)
    #fab.p[:,i,k]/=norm(fab.p[:,i,k])
    fab.str[i]=0

    end

  rs=fab.r[:,k];nbrs=fab.nbrs[:,:,k];grmob=fab.grmob;ngr=fab.ngr
  areas=fab.areas[:,:,k];p=fab.p[:,:,k]
  
 # fab.r[:,i]=advanceRadii(fab.r[:,i],fab.nbrs[:,:,i],fab.grmob,pars.dt,fab.epsdot[:,:,i],fab.ngr,fab.areas,fab.p[:,:,i])
#  rs_new=zeros(fab.ngr)
  vol=4/3.*pi.*rs.^3
  for i=2:ngr
        #partition
    for j in (1:i-1)[nbrs[1:i-1,i]]
      #get volume swept out be each boundary
      #this is relative to r[i]
      dVol=areas[i,j]*dt.*nggVelocity(rs[i],rs[j],grmob)
      dVol+=areas[i,j]*dt.*strEnVelocity(p[:,i],p[:,j],fab.epsdot[:,:,k]) ####!!!!!!!!!!!!!!!!!!!!!!!!
      dVol+=100*areas[i,j]*dt.*(fab.str[i,k]-fab.str[j,k])
      #print(dVol)
      vol[i]=vol[i]-dVol
      vol[j]=vol[j]+dVol
      if vol[j]<0
        vol[i]+=vol[j]
        vol[j]=0
        end
      if vol[i]<0
        vol[j]+=vol[i]
        vol[i]=0
        end
      end #j
    if rs[i]<0.4||fab.str[i]>0.02#r_crit  FIX !!!!!!!!!!!!!!!!!!!!!!
      T=0
      if rand()<prNucleation(T)#0.5#pr_nucleation
        nucleateGrain!(fab,i,k,vol)
        end
      end
    end #i
  fab.str+=dt  
  return (vol.*0.75/pi).^(1/3)

  end

function prNucleation(T)
  return 0.0#1
  end

##########################################      
function advanceRadius(this,rs,grmob,dt)
  if (this <= 0) | isnan(this)
    return 0.0 
  else
    dr=nggRate(this,rs,grmob,dt)
    new=this+dr
    if ((new<0) | (isnan(new)))
      return 0.0
    else
      return new
      end
    end
  end

function initAreas(rs,nbrs)
  ngr=length(rs)
  areas=zeros(ngr,ngr)
  for i=1:ngr
    areas[nbrs[:,i],i]=propAreas(rs[nbrs[:,i]])
    end
  areas=(areas+areas')/2
  return areas
  end

#get the velocity for ngg between two grains
#Important: OUTWARD from r1, don't fuck that up.
function nggVelocity(r1,r2,grmob)
  if ((r1>0) & (r2>0))
    mc=(1./r2-1./r1)./2
    return mc*grmob
  else return 0
    end
  end
#returns radius dt
function nggRate(this,rs,grmob,dt)
  pa=propAreas(rs)
  v=Array(Float64,length(rs))
  for i=1:length(rs)
    v[i]=nggVelocity(this,rs[i],grmob)
    end
  dV=sum(2/3*pi*(this*this*this-pa.*(this-v*dt).^3))
  if dV<0
    return -((abs(dV))^3)
    else
      return dV^(1/3)
    end
  end


function setConjugateNbrs(this::Int,nbr::Array{Float64,1},nbrs,nn,repl::Bool=false)
  for i=1:nn
    if nbr[i]==0
      nbr[i]=this
      return nbr
      end
    end
  i=Utils.randir(1,nn)
  nbrs[i]=this
  end

function siteRandomNbrs!(nbrs::Array{Int,2})
  nn,ngr=size(nbrs)
  for i=1:nn
    for j=1:ngr
        for k=1:nn
          if nbrs[i,j]==0
            di=Utils.diffrandi(j,1,ngr)
            if nbrs[k,di]==0
              nbrs[k,di]=j
              break
              end
            nbrs[1,di]=j
          end
        end
      end
    end
  end

makeRandomNbrs(ns::Int,ngr::Int,nn::Int)=makeRandomNbrs(zeros(Int64,nn,ns*ngr),ns,ngr,nn)

function makeSymNbrs(ns::Int,ngr::Int,pr)
  nbrM=fill(false,ngr,ngr,ns)
  for i=1:ns
    #change this#######################################################################
    for j=2:ngr
      for k=1:j-1  
        p=rand()
        if p<pr
          nbrM[j,k,i]=true
          nbrM[k,j,i]=true
          end
        end
      if all(x->x==false,nbrM[:,j,i])
          w=diffrandi(j,1,j-1)
          nbrM[j,w,i]=true
          nbrM[w,j,i]=true
          end
      end
    end
  return nbrM
  end



type GlobalPars{T<:Number,I<:Int}
  dt::T #timestep between velocity timesteps
  nrk::I #Number of timesteps to be taken per dt for Jeffery's eqn by RK4
  hrk::T #better be dt/nrk
  f::Function #Jeffery's eqn to supply to rk3
  function GlobalPars(dt,nrk,f)
    return new(dt,nrk,dt/nrk,f)
    end
  end


#function fabEv!(pars::GlobalPars,fab::FabricNGG,f)
#  for i=1:fab.ns
#    print("282")
#    fab.p[:,:,i]=nRK4(f,fab.ngr,fab.h,pars.nrk,pars.dt,fab.p[:,:,i],
#              fab.vort[:,:,i],fab.epsdot[:,:,i])
#    end
# for i=1:fab.ngr*fab.ns-1
#    print("287")
#    fab.r[i]=advanceRadius(fab.r[i],fab.r[filter(y->y!=0,fab.nbrs[:,i])],fab.grmob,pars.dt)
#    print("289\n",i)
#    end
#  return fab.r
#  end


#this is the closure that returns fabEvolve!,
#which evolves the fabric based on the timestep.
function fabricHelper(pars::GlobalPars,fab::FabricNGG,f::Function) #changed fab::AbstractFabric
  
  #function fabEvolve!(pars::GlobalPars,fab::Fabric3,f)
  function fabEvolve!(pars::GlobalPars,fab::FabricNGG,f)
    for i=1:fab.ns
      fab.p[:,:,i]=evFab!(
      #fab.r[:,i]=advanceRadii(fab.r[:,i],fab.nbrs[:,:,i],fab.grmob,pars.dt,fab.epsdot[:,:,i],fab.ngr,fab.areas,fab.p[:,:,i])
      fab.r[:,i]=advanceRadii(fab,i,pars.dt)
      end
#   for i=1:fab.ngr*fab.ns
#      fab.r[i]=advanceRadius(fab.r[i],fab.r[filter(y->y!=0,fab.nbrs[:,i])],fab.grmob,pars.dt)
#      end
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

function strEnVelocity(p1,p2,sigma)
  fac=1#10000
  e1=localSigmaEff([p1,p2],sigma,2)
  velocity=fac*(e1[2]-e1[1])
  return velocity
  end
#finds the effective stress on grains at one site.
#In addition, finds the local stress tensors for each grain
#at a site. 
function localSigmaEff(p,sigma,ngr)
  G=zeros(3,3)
  #first find local geometric tensors
  g=Array(Float64,3,3,ngr)
  m=Array(Float64,3,ngr)
  n=Array(Float64,3,ngr)
  for i=1:ngr
    g[9*i-8:9*i],m[3*i+1:3*i],n[3*i+1:3*i=localGeomTensor(p[3*i-2:3*i],sigma)
    G[1:9]+=g[9*i-8:9*i]
    end
  G=G/ngr
  sigmaE=zeros(ngr)
  for i=1:ngr
    #treat zeros!
    g[9*i-8:9*i]=g[9*i-8:9*i]./G[1:9].*sigma[1:9] #g= local sigma now
    
    sigmaE[i]=sqrt(abs(-1/3*secondInv(g[:,:,ngr]))) #be sure to convert
    #get stress tensor
    sigma[i

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
  T=sigma*p #read: T
  n=cross(p,T) #read: n
  n=n/norm(n) #read: n
  m=cross(n,p)
  lg=(m/norm(m))*p'
  return lg,m,n
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

function getRandc(n)
  c=rand(3,n)-0.5
  for i=1:n
     c[:,i]=c[:,i]/norm(c[:,i])
     end
  return c
  end
getRandc()=getRandc(100)

function evGroup!(p,sigma,dt)
  cpMat(v)=[0 v[3] -v[2];-v[3] 0 v[1]; v[2] -v[1] 0]
  function getRotM(x)
    (norm(x-[0,0,1])<1e-6)?(return eye(3)):pass
    v=cpMat(x)*z
    s=norm(v)
    c=x[3]
    V=cpMat(v)
    V[1,2]=V[1,2][1]
    V[1,3]=V[1,3][1]
    V[2,3]=V[2,3][1]
    V[2,1]=V[2,1][1]
    V[3,1]=V[3,1][1]
    V[3,2]=V[3,2][1]
     
    return eye(3)+V+(V*V)*(1-c)/s^2  
    end

  function getStrW(c)
    R=getRotM(c)
    sigmag=tensor2Voigt(R'*(voigt2Tensor(sigma))*R)
    Dgv[5]=sigmag[5];Dgv[4]=sigmag[4]
    Wg[2,3]=Dgv[4]
    Wg[3,2]=-Dgv[4]
    Wg[1,3]=Dgv[5]
    Wg[3,1]=-Dgv[5]
    Dg=voigt2Tensor(Dgv)
    Dg=R*Dg*R'
    Wg=R*Wg*R'
    print(Dg)
    return(Wg,Dg)
    end

  n=size(p)[2]
  W=zeros(3,3)
  D=zeros(3,3)
  Wg=zeros(3,3)
  Dgv=zeros(6)
  Dg=zeros(3,3)
  R=zeros(3,3)
  W=zeros(3,3)
  sigmag=zeros(3,3)
  z=[0,0,1]
  for i=1:n
    (Wg,Dg)=getStrW(p[:,i])
    W+=Wg;D+=Dg
    p[:,i]=grRotationRate(p[:,i],Wg)*dt+p[:,i]
    p[:,i]=p[:,i]/norm(p[:,i])
    Wg=zeros(3,3);Dg=zeros(3,3)
    end
  W/=n
  D/=n

  return (W,D)
  end #evgroup

#gets the rotation matrices
function getRotMHelper(p::Array{Float64,1},A::Array{Float64,2},
    A2::Array{Float64,2},R::Array{Float64,2})
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

end #module
