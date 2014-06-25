using Utils,jefferys,Plotting,PyCall
@pyimport matplotlib.pyplot as plt
ngr=300
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
vort_ss=zeros(3,3,ns)
vort_lg=zeros(3,3,ns)
vort_ss[2,3]=0.8
vort_ss[3,2]=-0.8
epsdot_ss=zeros(3,3,ns)
epsdot_lg=zeros(3,3,ns)
epsdot_lg[1,1,:]=-0.1#0.5
epsdot_lg[2,2,:]=0.2#-1
epsdot_lg[3,3,:]=-0.1#0.5
epsdot_ss[2,3,:]=0.8
epsdot_ss[3,2,:]=0.8

epsdot_lg_ss=epsdot_lg + epsdot_ss
vort_lg_ss = vort_ss

dt=1e-1
#5e-3
nrk=10
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




#fab=Fabric{Float64,Int64}(coors,p,ngr,ns,h,C,vort,epsdot)
#fabE=fabricHelper(pars,fab,jefferysRHS)
#fabE(pars,fab,jefferysRHS)
epsdot = epsdot_ss
vort = vort_ss

sv=Array(Float64,0)
for i=1:100
 fabE(pars,fab,jefferysRHS)
 x=svd(fab.p[:,:,1])[2]
 sv=append!(sv,[min(x)/norm(x)])
 end
schmidtPlot(fab.p);plt.show()



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




