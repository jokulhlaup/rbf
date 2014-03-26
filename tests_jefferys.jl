using Utils,jefferys,Plotting,PyCall
@pyimport matplotlib.pyplot as plt
ngr=30
ns=1
coors=rand(3,ns)
p=rand(3,ngr,ns)-0.5


function jefferysRHS(c,vort,epsdot,dt)
  return vort*dt*c+(-epsdot*c+(c'*epsdot*c)[1]*c)*dt
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
vort[3,2]=-0.4
vort[2,3]=0.4
epsdot=zeros(3,3,ns)
epsdot[1,1,:]=-0.05#0.5
epsdot[2,2,:]=0.1#-1
epsdot[3,3,:]=-0.05#0.5
epsdot[2,3,:]=0.4
epsdot[3,2,:]=0.4

dt=0.00005
#5e-3
nrk=10
f=jefferysRHS
pars=GlobalPars{Number,Int64}(dt,nrk,f)
#radius velocity length/time * dt in units length
#nggVelocity units 1/length*
nn=10
av_radius=1
r=rand(ngr*ns)
#fab=FabricNGG{Float64,Int64}(coors,p,ngr,ns,h,C,vort,epsdot,nn,nbrs,r)
fab=consFabricNGG(coors,p,ngr,ns,h,C,vort,epsdot,nn,av_radius)
fabE=fabricHelper(pars,fab,jefferysRHS)




#fab=Fabric{Float64,Int64}(coors,p,ngr,ns,h,C,vort,epsdot)
#fabE=fabricHelper(pars,fab,jefferysRHS)
#fabE(pars,fab,jefferysRHS)
sv=Array(Float64,0)
for i=1:100
 fabE(pars,fab,jefferysRHS)
 x=svd(fab.p[:,:,1])[2]
 sv=append!(sv,[min(x)/norm(x)])
 end




function obj_epsdot(epsd)
  



pl=schmidtPlot(fab.p)
plt.show()

####
#Some unit tests
###
rs=sort(rand(10))
assert(sum(jefferys.propAreas(rs))==1)
assert(nggVelocity(1,1,1)==0)
assert(nggRate(1,2,0.1,0.0001)<0)




