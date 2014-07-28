#test alignFabrics
p1=Utils.getRandc(10)
p1[3,:]=0.
for i=1:10
  p1[:,i]=p1[:,i]/norm(p1[:,i])
  end
p2=rotp(pi/2,p1)
(xmin,fmin,pmin)=alignFabrics(p1,p2)
@assert sum(abs(p1-pmin)/30)<0.01


