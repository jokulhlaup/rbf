
using jefferys

#Get 100 initial particles in (theta,phi)
#This converges to near SS at around
#t=100 with step size h=0.1

#n=theta->[sin(theta[1])*cos(theta[2]),sin(theta[1])*sin(theta[2]),cos(theta[1])]
#f=theta->((vort*n(theta))[2:3])-((epsdot*n(theta))[2:3]-(n(theta)'*epsdot*n(theta))[1]*n(theta)[2:3])

function f(theta)
  epsdot=eye(3)
  vort=zeros(3,3)
  epsdot[3,3]=-1
  epsdot[1,1]=100#sqrt(0.5)
  epsdot[2,2]=0.0
  n=[sin(theta[1])*cos(theta[2]),sin(theta[1])*sin(theta[2]),cos(theta[1])]
  return (((vort*n)[2:3])-((epsdot*n)[2:3]-((n')*epsdot*n)[1]*n[2:3]))
end
#using anonymous
n=100
theta=rand(2,n)
coors=rand(3)
#R=[cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)]
#R=[1 0 0;0 cos(theta) -sin(theta);0 sin(theta) cos(theta)]
 C=zeros(6,6)
 C[5,5]=1


#####################
#Tests for fabric
#####################
dt=0.1
nrk=10
hrk=dt/nrk
fab=FabricPt(theta,coors,n)
par=GlobalPars(dt,nrk,hrk,f)
jefferys.fabEvolve!(fab,par);

