module Wrap
@pyimport matplotlib.pyplot as plt 
ngr=100
ns=10
coors=rand(3,ns)
p=rand(3,ngr,ns)-0.5


function jefferysRHS(c,vort,epsdot,dt)
  return (epsdot*c-(c'*epsdot*c)[1]*c)*dt
    end 

#normalize!
for i=1:ngr
  for j=1:ns
     p[:,i,j]=p[:,i,j]/norm(p[:,i,j])
   end 
 end 
h=1.0
C=Array(Float64,6,6,ns)
vort=zeros(3,3,ns)
epsdot=zeros(3,3,ns)
epsdot[1,1,:]=1
epsdot[2,2,:]=1
epsdot[3,3,:]=-2


dt=5e-3
nrk=10
f=jefferysRHS
pars=GlobalPars{Float64,Int64}(dt,nrk,f)

nn=10
av_radius=1
r=rand(ngr*ns)
#fab=FabricNGG{Float64,Int64}(coors,p,ngr,ns,h,C,vort,epsdot,nn,nbrs,r)
fab=consFabricNGG(coors,p,ngr,ns,h,C,vort,epsdot,nn,av_radius)
fabEv!(pars,fab,jefferysRHS)
fabE=fabricHelper(pars,fab,jefferysRHS)
fab.p=fabE(pars,fab,jefferysRHS)

function objective(fab::FabricNGG,

