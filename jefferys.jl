module jefferys
using Utils,Distributions
export consFabricNGG,getRotM,Fabric,Fabric2,FabricNGG,genrFT,makeRandomNbrs!,fabEv!,advanceRadius,GlobalPars,AbstractFabric,solveJefferys,rk4,nRK4,rotC,jefferysRHS,fabricHelper,propAreas,getC,polygonize!,fp_soft_C
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
  C::Array{T,3} #viscosity matrix
  vort::Array{T,3} #vorticity
  epsdot::Array{T,3} #strain rate
  #Fabric(coors,p,ngr,h,ns,C=zeros)=new(coors,p,ngr,h,ns,C)
  #stencil::Array{T,1}
  nn::I
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
    return new(coors,p,ngr,h,ns,C,vort,epsdot,nn,nbrs,r)
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
  temp::Float64
  end))

function consFabricNGG(coors,p,ngr,ns,h,C,vort,epsdot,nn,av_radius,temp)
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
  size(epsdot)==(3,3,ns)?nothing:error("Dimension mismatch in 'epsdot'")
  return FabricNGG{Float64,Int64}(coors,p,ngr,h,ns,C,vort,epsdot,nn,nbrs,areas,r,grmob,pr_nuc,nuc_vol,str,temp)
  end

######################
############
#####################
##########################################
macro nanch(test_var)
  quote
  if any(isnan,$test_var)
    global ERR_VAR=$test_var
    error("NaN in ", $test_var)
  end
  end
  end

macro nanch2(test_var,name)
  quote
    if any(isnan($test_var))
      error("NaN in ", $test_var," ",$name)
      end
    end
  end

velgrad(sigma,Sx,Sy)=Sx*ddot(Sx,sigma)+Sy*ddot(Sy,sigma)#+Sz*ddot(Sz,sigma)
ddot(A,B)=trace(A*B')

function get_rss_softness(fab,sigma,stress_fac,k)
  n=size(fab.p,2)
  rss_0=Array(Float64,n)
  softness=Array(Float64,n)
  rst=Array(Float64,6,n)
  #the resolved strain tensor (voigt)
  for i=1:n
    R=getRotM(fab.p[:,i,k])
    rst[:,i]=tensor2Voigt(R*sigma*R')
    rss_0[i]=sqrt(0.5*(rst[4,i]^2+rst[5,i]^2)) 
    end
  for i=1:n
    softness[i]=((1-stress_fac) + dot(stress_fac*fab.areas[:,i],rss_0)/rss_0[i]) 
    end
  return (rst,rss_0,softness)
  end

#get the compatible bulk sigma and softness
function fp_soft_C(fab,stress_fac,k)
  n=size(fab.p,2)
  C=getC(fab,k)
  sigma=voigt2Tensor(getC(fab,1)*tensor2Voigt(fab.epsdot[:,:,k]))
  (rst,rss_0,softness)=get_rss_softness(fab,sigma,stress_fac,k);
  misfit=1.
  while (misfit > 1e-6)
    sigma_old=sigma
    (rst,rss_0,softness)=get_rss_softness(fab,sigma,stress_fac,k)      
    sigma=voigt2Tensor(getC(fab,softness,1)*tensor2Voigt(fab.epsdot[:,:,k]))
    misfit=sum(abs((sigma-sigma_old)/Utils.secondInv(sigma)))
    end
  return (rst,rss_0,softness,sigma)
  end

function thorRot!(fab,pars,k,dt,stress_fac)
  #get bulk sigma
  #bulk sigma in voigt; 
  bulk_sigma=voigt2Tensor(getC(fab,1)*tensor2Voigt(fab.epsdot[:,:,k]))
  # bulk_sigma=zeros(3,3);
  #bulk_sigma[1,1]=1
  #bulk_sigma[2,2]=1
  #bulk_sigma[3,3]=-2

  #now, get RSS on the basal plane for each grain
  
  n=size(fab.p,2)
  b0=zeros(3,n);
  b1=zeros(3,n);
  rss=Array(Float64,n)
  rss_0=Array(Float64,n)
  softness=Array(Float64,n)
  S0=Array(Float64,3,3,n);S1=Array(Float64,3,3,n)
  vort=zeros(3,3,n)
  rst=Array(Float64,6,n);
  sigma=voigt2Tensor(inv(getC(fab,k))*tensor2Voigt(fab.epsdot[:,:,k]))
  eff_stress=sqrt(0.5*abs(Utils.secondInv(sigma)))
  C=zeros(6,6)
  C[5,5]=1.
  C[6,6]=0.01
  C[1,1]=0.01
  C[2,2]=0.01
  C[3,3]=0.01
  C[4,4]=1.

  for i=1:n
    if i==50;print('C',C);end
    R=getRotM(fab.p[:,i,k])
    #the resolved strain tensor (voigt)
    rst[:,i]=C*tensor2Voigt(R*sigma*R')
    rss_0[i]=sqrt(0.5*(rst[4,i]^2+rst[5,i]^2)) 
#    vort[1,3,i]=rst[5,i]
#    vort[1,2,i]=rst[6,i]
#    vort[2,1,i]=-rst[6,i]
#    vort[2,3,i]=rst[4,i]
#    vort[3,1,i]=-rst[5,i]
#    vort[3,2,i]=-rst[4,i]
    #vort[:,:,i]=R'*vort[:,:,i]*R
#    b0[3,i]=0;
#    b0[1,i]=1;
#    b0[2,i]=-b0[1,i]*fab.p[1,i,k]/fab.p[2,i,k];
#    b0[:,i]=b0[:,i]/norm(b0[:,i])
#    b1[:,i]=cross(b0[:,i],fab.p[:,i,k])
    #S0[:,:,i]=(b0[:,i]*fab.p[:,i,k]' )'
    #S1[:,:,i]=(b1[:,i]*fab.p[:,i,k]')'
    #rss_0[i]=norm(ddot(S0[:,:,i],sigma)*b0[:,i] + ddot(S1[:,:,i],sigma)*b1[:,i])
    #traction=sigma*fab.p[:,i,k]
    #norm_str=traction*p
    #rss=sqrt(dot(traction,traction)-dot(norm_str,norm_str))
    end
  for i=1:n
    softness[i]=((1-stress_fac) + dot(stress_fac*fab.areas[:,i],rss_0)/rss_0[i])
    #gamma_0=ddot(S0[:,:,i],sigma)*softness[i]
    #gamma_1=ddot(S1[:,:,i],sigma)*softness[i]
#    v_grad=S0[:,:,i]*gamma_0+S1[:,:,i]*gamma_1
    #fab.p[:,i,k]+=(v_grad-v_grad')*fab.p[:,i,k]*dt
    if(i==50);print("delta, ", vort[:,:,i]*fab.p[:,i,k]*dt);end
    #fab.p[:,i,k]-=vort[:,:,i]*fab.p[:,i,k]*dt

    #fab.p[:,i,k]=-softness[i]*rk4(pars.f,fab.ngr,fab.p[:,i,k],vort[:,:,i],voigt2Tensor(rst[:,i]),dt);
    rotf=Utils.fisher_rot_mat(100.)
    rvort=rotf'*fab.vort[:,:,k]*rotf
    repsdot=rotf'*fab.epsdot[:,:,k]*rotf
    dx=rk4(fab.ngr,fab.p[:,i,k],rvort,repsdot,dt,softness[i]);
    fab.p[:,i,k]+=dx
    poly_stress_ratio=0.09;#0.065;
    poly_angle_diff=0.087266
    if abs((rss_0[i]/eff_stress)<poly_stress_ratio)
        polygonize!(fab,dx,i,k,poly_angle_diff)
    end

#    fab.p[:,i,k]+=rand(Distributions.Gaussian(0,0.05),3)
    fab.p[:,i,k]/=norm(fab.p[:,i,k])
  end
end
#function viscRot!(fab,k,dt,stress_fac)
function polygonize!(fab,dx,i,k,poly_angle_diff)
   for l=1:length(fab.r)
       if fab.r[l,k]<1e-8
          dxn=-dx/norm(dx)
          fab.p[:,l,k]=Utils.rotate_towards!(fab.p[:,i,k],dxn,-poly_angle_diff)
          fab.p[:,l,k]+=rand(Distributions.Gaussian(0,0.15),3)
          fab.p[:,l,k]/=norm(fab.p[:,l,k])
          fab.r[l,k]=(fab.r[i,k]^3*0.5)^(1/3)
          fab.r[i,k]=fab.r[l,k]
          print("polygon!")
          return true;
       end
   end
   return false
end

#function advanceRadii(rs,nbrs,grmob,dt,sigma,ngr,areas,p,pr_nucleation,nuc_vol)
function advanceRadii(fab::FabricNGG,k,dt)
  
  #nucleates a new grain at position (i,j,k)
  function nucleateGrain!(fab,i,k,vol)
    jd=binBoolInd(fab.areas[:,i,k],>,ngr)
    vol[jd]-=fab.nuc_vol
    vol[i]=fab.nuc_vol

    @nanch(vol)
    if vol[jd]<0
      vol[i]+=vol[jd]
      vol[jd]=0.
      end
    @nanch(vol)
    fab.p[:,i,k]=getRandOrient()
    fab.str[i,k]=0.

    end
  function getRandOrient()     
    azimuth=rand()*2.0*pi  
    zenith=rand()*pi/2.0
    return [sin(zenith)*cos(azimuth),sin(zenith)*sin(azimuth),cos(zenith)]#rand(3)
    #fab.p[:,i,k]/=norm(fab.p[:,i,k])
    end

  rs=fab.r[:,k];nbrs=fab.nbrs[:,:,k];grmob=fab.grmob;ngr=fab.ngr
  areas=fab.areas[:,:,k];p=fab.p[:,:,k]
  rss_0=Array(Float64,ngr)
  
  C=0.01*eye(6);C[5,5]=1.;C[4,4]=1
  sigma=voigt2Tensor(inv(getC(fab,k))*tensor2Voigt(fab.epsdot[:,:,k]))
  eff_stress=sqrt(0.5*abs(Utils.secondInv(sigma)))
  
  for i=1:ngr
    R=getRotM(fab.p[:,i,k])
    #the resolved strain tensor (voigt)
    rst=C*tensor2Voigt(R*sigma*R')
    rss_0[i]=sqrt(0.5*(rst[4]^2+rst[5]^2)) 
#    fab.p[:,i,k]=polygonize!(fab,indmax_eigval,eigenstoff,eff_stress,rss_0[i],i,k)
    end

  @nanch(rs)
 # fab.r[:,i]=advanceRadii(fab.r[:,i],fab.nbrs[:,:,i],fab.grmob,pars.dt,fab.epsdot[:,:,i],fab.ngr,fab.areas,fab.p[:,:,i])
#  rs_new=zeros(fab.ngr)
  vol=4/3.*pi.*rs.^3
  for i=1:ngr
    if any(isnan,p[:,i])
      p[:,i]=getRandOrient()
      end
    end
  if any(isnan,sigma)
      sigma=zeros(size(sigma))
  end
  @nanch(sigma)
  for i=2:ngr
        #partition
    if any(isnan,p[:,i])
      p[:,i]=getRandOrient()
      end
    for j in (1:i-1)[nbrs[1:i-1,i]]
      #get volume swept out be each boundary
      #this is relative to r[i]
      dVol=-areas[i,j]*dt.*nggVelocity(rs[i],rs[j],grmob)
      @nanch(dVol)
      if any(isnan,p[:,j])
        p[:,j]=getRandOrient()
        end
      @nanch(p[:,j])
      @nanch2(p[:,i],"p[:,i]")
      @nanch2(sigma,"sigma")
      @nanch2(fab.str,"string")
      #########
      ##
      dVol+=areas[i,j]*dt.*grmob*((rss_0[i]-rss_0[j]) + 
          strEnVelocity(p[:,i],p[:,j],sigma[:,:,k],fab.str[i,k],fab.str[j,k])) 
      ####!!!!!!!!!!!!!!!!!!!!!!!!
      @nanch2(dVol,"dVol")
#      dVol+=100*areas[i,j]*dt.*(fab.str[i,k]-fab.str[j,k])
      #print(dVol)
      vol[i]=vol[i]-dVol
      vol[j]=vol[j]+dVol
      @nanch2(vol,"vol")
      (isnan(rs[i])||isnan(rs[j]))?error("isNaN: ", i,", ", j):nothing
      if vol[j]<0
        vol[i]+=vol[j]
        vol[j]=0.
        end
      if vol[i]<0
        vol[j]+=vol[i]
        vol[i]=0.
        end
      end #j
    if rs[i]<1e-4#r_crit
      if rand()<prNucleation(fab.temp,dt)#0.5#pr_nucleation
        nucleateGrain!(fab,i,k,vol)
        println("nucleated!!!!!!!!!!")
        end
      end
    end #i
  
  fab.str+=dt*sqrt(abs(secondInv(fab.epsdot[:,:,k])))  
  return (vol.*0.75/pi).^(1/3)

  end

function prNucleation(A,b,T)
  return A*exp(b*(T))#1
  end

prNucleation(T,dt)=prNucleation(1.,0.03/dt,T)
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

#Modification of ODE4 from package ODE
function nRK4(f,ntimes::Int,h::Number,m::Number,dt,p,vort,epsdot)
  for i=1:ntimes
     p[:,i]=jefferys.rk4(h,m,p[:,i],vort[:,:],epsdot[:,:],dt,m)
     end
  return p
  end

function rk4(f::Function,n::Int64,x,vort_b,epsdot_b,dt)
   #vort=deepcopy(epsdot); vort[3,1:2]=0
   #vort[2,1]=0; vort=vort-vort'
   for i=1:n
      k1=f(x,vort,epsdot,dt)
      k2=f(x+k1*dt/2,vort,epsdot,dt)
      k3=f(x+k2*dt/2,vort,epsdot,dt)
      k4=f(x+k3*dt/2,vort,epsdot,dt)
      x+=(1/6)*dt*(k1+2*k2+2*k3+k4)
      end
   return x
   end


function ejefferys(x,vort,epsdot,dt)
  res=zeros(3)
  for j=1:3
    for i=1:3
      res[i]-=epsdot[i,j]*x[j]
      end
    end
  epsc2=res[1]*x[1]+res[2]*x[2]+res[3]*x[3]
  for j=1:3
    @simd for i=1:3
      @inbounds res[i]+=vort[i,j]*x[j]
     end
    end
  res=(res-epsc2*x)*dt
  end
   
#function ejefferys(x,vort,epsdot,dt)
#   res=zeros(3)
#   for j=1:3
#       @simd for i=1:3
#         @inbounds res[i]=-epsdot[i,j]*x[j]
#      end
#   end
#   epsc2=res[1]*x[1]+res[2]*x[2]+res[3]*x[3]
#   for j=1:3
#      @simd for i=1:3
#        @inbounds res[i]+=vort[i,j]*x[j]
#      end
#   end
   
   #res2=zeros(3)
   #res2[1]=epsdot[1,3]*x[3]
   #res2[2]=epsdot[2,3]*x[3]
   #res2[
   ##res[3]-=epsdot[2,3]*x[2]
   #res[3]-=epsdot[1,3]*x[1]
   #res=(res-epsc2*x)*dt
   #vt=zeros(3,3); vt[1,2:3]=epsdot[1,2:3]
   #vt[2,3]=epsdot[2,3];vt-=vt'
   #res=0.8*res+0.2*(vt*x-epsdot*x+
#   return res - epsc2*x
#end
       
#function ejefferys(x,vort,epsdot,dt)
#   res=zeros(3)
#   for j=1:3
#      for i=1:3
#         res[i]-=epsdot[i,j]*x[j]
#      end
#   end
#   epsc2=-res[1]*x[1]-res[2]*x[2]-res[3]*x[3]
#   for j=1:3
#      @simd for i=1:3
#        @inbounds res+=vort[i,j]*x[j]
#      end
#   end
#   res=(res+epsc2*x)*dt
#end
 
function rk4(n::Int64,x,vort,epsdot_b,dt,softness)
  a=(softness/6)*dt
  dx=0.
  rot_mat=Utils.rotate_mat_rodriguez(x)
  epsdot=epsdot_b
  epsdot=rot_mat*epsdot_b*rot_mat'
  epsdot13=epsdot[1,3];epsdot23=epsdot[2,3]
  epsdot=zeros(3,3);epsdot[1,3]=epsdot13
  epsdot[2,3]=epsdot23
  epsdot+=epsdot'
  epsdot=rot_mat'*epsdot*rot_mat
#
  for i=1:n
    k1=ejefferys(x,vort,epsdot,dt)
    k2=ejefferys(x+k1*dt/2,vort,epsdot,dt)
    k3=ejefferys(x+k2*dt/2,vort,epsdot,dt)
    k4=ejefferys(x+k3*dt/2,vort,epsdot,dt)
    dx+=(softness/6)*dt*(k1+2*k2+2*k3+k4)
    end
  return dx
  end

#2function rk4(n::Int64,x,vort,epsdot_b,dt,softness)
#   a=(softness/6)*dt
#   dx=zeros(3)
#   #xr=[0,0,1]
##   rot_mat=Utils.rotate_mat_rodriguez(x)
#   epsdot=epsdot_b
##   epsdot=rot_mat*epsdot_b*rot_mat'
##   epsdot13=epsdot[1,3];epsdot23=epsdot[2,3]
##   epsdot=zeros(3,3);epsdot[1,3]=epsdot13
##   epsdot[2,3]=epsdot23
##   epsdot+=epsdot'
##   epsdot=rot_mat'*epsdot*rot_mat
#   for i=1:n
#      k1=ejefferys(x,vort,epsdot,dt)
#      k2=ejefferys(x+k1*dt/2,vort,epsdot,dt)
#      k3=ejefferys(x+k2*dt/2,vort,epsdot,dt)
#      k4=ejefferys(x+k3*dt/2,vort,epsdot,dt)
#      dx+=(softness/6)*dt*(k1+2*k2+2*k3+k4)
#      end
##   dx=rot_mat*dx   
#   return dx
#   end


function jefferysRHS(c,vort,epsdot,dt)
#  return (vort*c + epsdot*c-(c'*epsdot*c)[1]*c)*dt
   return (-vort*c + epsdot*c -(c'*epsdot*c)[1]*c)*dt
  end
#Replace this so its rotC(R)
function rotC(R)
  C = zeros(6,6)
  C[5,5]=1 
  C[4,4]=1
  rotC(R,C)
  end
 
function rotC(R,C)
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
         R[1,2].*R[2,3]+R[1,3].*R[2,2] R[1,3].*R[2,1]+ R[1,1].*R[2,3] R[1,1].*R[2,2]+R[1,2].*R[2,1]] ; 
  K = [ K1  2*K2 ; 
        K3   K4   ] ;
  C = K * C * transpose(K) 
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
  #this is the function that actually does the rotation.  
  function fabEvolve!(pars::GlobalPars,fab::Fabric,f) #jefferys equation
    for i=1:fab.ns*fab.ngr
      if (fab.r[i] > 0)
        fab.p[:,i]=rk4(f,fab.ngr,fab.h,pars.nrk,pars.dt,fab.p[:,i,:],
          fab.vort[:,:,i],fab.epsdot[:,:,i])
        end
      end
      return fab.p
    end
  
  #function fabEvolve!(pars::GlobalPars,fab::Fabric3,f)
  function fabEvolve!(pars::GlobalPars,fab::FabricNGG,f)
    for i=1:fab.ns
      #update the orientations
      thorRot!(fab,pars,i,pars.dt,1.)
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

function strEnVelocity(p1,p2,sigma,str1,str2)
  if all(sigma.==0)
      return 0.
  end
  fac=1.
  e1=localSigmaEff([p1,p2],sigma,2)
  @nanch(e1)
  velocity=fac*(e1[2]-e1[1])+(str2-str1)
  @nanch(velocity)
  return velocity
  end
#finds the effective stress on grains at one site.
#In addition, finds the local stress tensors for each grain
#at a site. 
function localSigmaEff(p,sigma,ngr)
  if all(sigma.==0)
      return 0.
  end
  G=zeros(3,3)
  @nanch2(sigma,"sigma492")
  @nanch2(p,"p 493")
  if all(x->x==0,sigma)
      error("sigma ekwals zero")
  end
  #first find local geometric tensors
  g=zeros(Float64,3,3,ngr)
  m=zeros(Float64,3,ngr)
  n=zeros(Float64,3,ngr)
  for i=1:ngr
    o=localGeomTensor(p[3*i-2:3*i],sigma)
    g[9*i-8:9*i]=o[1]
    m[3*i-2:3*i]=o[2]
    n[3*i-2:3*i]=o[3]
    G[1:9]+=g[9*i-8:9*i]
    end
  G=G/ngr
  @nanch2(G,"G")
  @nanch2(g,"g1")
  sigmaE=zeros(ngr)
  for i=1:ngr
    #treat zeros!
    g[9*i-8:9*i]=g[9*i-8:9*i]./G[1:9].*sigma[1:9] #g= local sigma now
    @nanch2(g,"g")
    
    sigmaE[i]=sqrt(abs(-1/3*secondInv(g[:,:,ngr]))) #be sure to convert
    @nanch2(sigmaE,"sigmae")
    #get stress tensor
    #sigma[i

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
  return (lg,m,n)
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
function getRotM(x)
  (norm(x-[0,0,1])<1e-6)?(return eye(3)):nothing
  cpMat(v)=[0 v[3] -v[2];-v[3] 0 v[1]; v[2] -v[1] 0]
  z=[0,0,1]
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

function getC(fab,softness,k)
  tot_vol=sum(fab.r[:,:,k].^3) #unscaled b/c it is divided anyway
  #change this such that accepts arbitrary C_p
  C=zeros(6,6)
  C_p=100.*eye(6)
  C_p[4,4]=1.;C_p[5,5]=1.
  for i=1:fab.ns*fab.ngr
    C_p[4,4]=1.0*softness[i];C_p[5,5]=1.0*softness[i];
    R=getRotM(fab.p[3*i-2:3*i])
    C+=rotC(R,C_p).*(fab.r[i]^3)
    end
  C/=tot_vol
  end

function getC(fab::AbstractFabric,k)
  tot_vol=sum(fab.r[:,:,k].^3) #unscaled b/c it is divided anyway
  #change this such that accepts arbitrary C_p
  C=zeros(6,6)
  C_p=100.*eye(6)
  C_p[4,4]=1.;C_p[5,5]=1.
  for i=1:fab.ns*fab.ngr
    R=getRotM(fab.p[3*i-2:3*i])
    C+=rotC(R,C_p).*(fab.r[i]^3)
    end
  C/=tot_vol
  end

#Main driver routine to get the viscosity.
#refactor this by 
#at step zero, generate a closure equiv to fabEvolve! +params
#then at each timestep, mutate the closure. (can you do that
#without pushing a new copy onto the stack?)
end #module
