using jefferys,Plotting,PyCall
@pyimport matplotlib.pyplot as plt
ngr=100
ns=10
coors=rand(3,ns)
p=rand(3,ns,ngr)-0.5
#normalize!
for i=1:ns
  for j=1:ngr
     p[:,i,j]=p[:,i,j]/norm(p[:,i,j])
   end
 end
h=1
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

fab=Fabric{Float64,Int64}(coors,p,ngr,ns,h,C,vort,epsdot)
fabE=fabricHelper(pars,fab,jefferysRHS)
fabE(pars,fab,jefferysRHS)
for i=1:10
 fab.p=fabE(pars,fab,jefferysRHS)
 end

pl=schmidtPlot(fab.p)
plt.show()


