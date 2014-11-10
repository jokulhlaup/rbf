using jefferys

#Get 100 initial particles in (theta,phi)
#This converges to near SS at around
#t=100 with step size h=0.1

#n=theta->[sin(theta[1])*cos(theta[2]),sin(theta[1])*sin(theta[2]),cos(theta[1])]
#f=theta->((vort*n(theta))[2:3])-((epsdot*n(theta))[2:3]-(n(theta)'*epsdot*n(theta))[1]*n(theta)[2:3])
epsdot1=eye(3)
vort1=zeros(3,3)
epsdot1[3,3]=0.
epsdot1[1,1]=-1.#sqrt(0.5)
epsdot1[2,2]=0.

const epsdot=epsdot1
const vort=vort1

function f(p)
  m=(vort*p + (epsdot*p)-((p')*epsdot*p)[1]*p)#+(vort*p)
  return m
  end
#using anonymous
n=100
p=rand(3,n)
for i=1:n
  p[:,i]=p[:,i]/norm(p[:,i])
  end
coors=rand(3)
#R=[cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)]
#R=[1 0 0;0 cos(theta) -sin(theta);0 sin(theta) cos(theta)]
 C=zeros(6,6)
 C[5,5]=1


#####################
#Tests for fabric
#####################
dt=0.01
nrk=10
hrk=dt/nrk
fab=FabricPt(p,coors,n)
par=GlobalPars(dt,nrk,hrk,f)
jefferys.fabEvolve!(fab,par);

