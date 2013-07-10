
using jefferys

#Get 100 initial particles in (theta,phi)
#This converges to near SS at around
#t=100 with step size h=0.1
theta=rand(2,100000)
theta[1,:]=theta[1,:]*pi
theta[2,:]=(theta[2,:]-0.5)*2.0*pi
epsdot=eye(3)
vort=zeros(3,3)
epsdot[3,3]=-1
epsdot[1,1]=100#sqrt(0.5)
epsdot[2,2]=0.0
a=rand(size(theta[:,1]))

n=theta->[sin(theta[1])*cos(theta[2]),sin(theta[1])*sin(theta[2]),cos(theta[1])]
f=theta->((vort*n(theta))[2:3])-((epsdot*n(theta))[2:3]-(n(theta)'*epsdot*n(theta))[1]*n(theta)[2:3])

#using anonymous 

