module Constructors

using Utils,Plotting,PyCall,jefferys
export jefferysRHS,mkFab
@pyimport matplotlib.pyplot as plt

function jefferysRHS(c,vort,epsdot,dt)
 return vort*dt*c+(-epsdot*c+(c'*epsdot*c)[1]*c)*dt
   end

function mkFab(ngr)
ns=1
coors=rand(3,ns)
p=rand(3,ngr,ns)-0.5


#normalize!
for i=1:ngr
  for j=1:ns
     p[:,i,j]=p[:,i,j]/norm(p[:,i,j])
   end
 end
h=1.0
C=Array(Float64,6,6,ns)
vort=zeros(3,3,ns)
vort[2,3]=0.8
vort[3,2]=-0.8

epsdot=zeros(3,3,ns)
#epsdot[1,1,:]=-0.1#0.5
#epsdot[2,2,:]=0.2#-1
#epsdot[3,3,:]=-0.1#0.5
#epsdot[2,3,:]=0.8
#epsdot[3,2,:]=0.8

dt=1e-1
#5e-3
nrk=1000
f=jefferysRHS
pars=GlobalPars{Number,Int64}(dt,nrk,f)
#radius velocity length/time * dt in units length
#nggVelocity units 1/length*
nn=10
av_radius=1
temp=-100.
r=rand(ngr*ns)
#fab=FabricNGG{Float64,Int64}(coors,p,ngr,ns,h,C,vort,epsdot,nn,nbrs,r)
fab=consFabricNGG(coors,p,ngr,ns,h,C,vort,epsdot,nn,av_radius,temp)
fabE=fabricHelper(pars,fab,jefferysRHS)
return (fabE,fab,pars)
end

mkFab()=mkFab(30)
end #module

