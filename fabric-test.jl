
using jefferys

#Get 100 initial particles in (theta,phi)
theta=rand(2,100)
theta[:,1]=theta[:,1]*pi
theta[:,2]=(theta[:,2]-0.5)*2.0*pi

a=rand(size(theta[:,1]))

n=theta->[sin(theta[1])*cos(theta[2]),sin(theta[1])*sin(theta[2]),cos(theta[1])]
f=theta->((vort*n(theta))[2:3])-((epsdot*n(theta))[2:3]-(n(theta)'*epsdot*n(theta))[1]*n(theta)[2:3])

function nRK4(f,ntimes,h,m(
  @parallel (+) for i=1:ntimes
     jefferys.solveJefferys(eye(3),eye(3),theta[:,rand[1:ntimes],0.1,100)
     end
