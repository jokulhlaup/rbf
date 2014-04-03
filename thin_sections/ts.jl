using Loess,Grid,PyCall
@pyimport matplotlib.pyplot as plt

function wrapper(ages,par,ts_svs,fab,pars,jefferysRHS)
  t_c=p[1]
  fab.epsdot[1,1]=par[2]
  fab.epsdot[2,2]=par[3]
  fab.epsdot[3,3]=par[4]
  fab.epsdot[2,3]=par[5]
  fab.epsdot[1,3]=par[6]
  fab.epsdot[1,2]=par[7]
  symmetrize!(fab.epsdot)
  fab.vort[3,2]=-par[8]
  fab.vort[2,3]=par[8]
  fab.epsdot[2,3]=par[8]
  fab.epsdot[3,2]=par[8]
  function objective(p)
    sv=Array(Float64,0)
    dages=ages[2:end]-ages[1:end-1]
    for i=1:length(dages)
      pars.dt=t_c*dage
      fabE(pars,fab,jefferysRhs)
      sv[:,i]==svd(fab.p[:,:,1])[2]
      end
    sse=sum((ts_svs-sv).^2)
    return sse
    end
  return objective
  end
  
function readif(rf)
  c=Array(Float64,0)
  for i=1:length(rf)
    try  append!(c,float64(split(replace(rf[i],r"[\t\r]","\r"))[2:3]))  #append!(c,float64(split(rf[i],'\t')[2:3]))
      catch
        continue
      end
    end
  return reshape(c,(2,int(length(c)/2)))
  end

cd("thin_sections/C-axisdatabase")
dr=readdir()
c=Dict()
pd=Dict()
svs=Array(Float64,length(dr))
x=Array(Float64,3,length(dr))
for d=1:length(dr)
  cd(dr[d])
  rf=readdlm("c-axes.txt",'\n')
  i=int(dr[d])
  c[i]=readif(rf)
  #convert to rads
  #c[i][2,:] is theta (azimuth angle)
  c[i]=c[i]*pi/180
  #get regular pi
  #c[i][1,:]=acos(sin(c[i][2,:]).*sin(c[i][1,:]))
  pd[i]=Array(Float64,3,length(c[i][1,:]))
  pd[i][1,:]=cos(c[i][1,:])
  pd[i][2,:]=sin(c[i][1,:]).*cos(c[i][2,:])
  pd[i][3,:]=sin(c[i][1,:]).*sin(c[i][2,:])
  for j = 1:length(pd[i])/3
    if pd[i][3,j] < 0
    pd[i][1,j] = - pd[i][1,j]  ;
    pd[i][2,j] = -pd[i][2,j];
    pd[i][3,j] = -pd[i][2,j] ;
      end
    end
  x[:,d]=svd(pd[i])[2]
  svs[d]=minimum(x[:,d])/norm(x[:,d])
  cd("..")
  x[:,d]=sort(x[:,d])/norm(x[:,d])
  end

cd("..")
cs=readdlm("timescale.csv",',')



#now interpolate depth-age
dr=sort(float(dr))

depth_temps=readdlm("temp.TXT")
predicted_temps=loess(depth_temps[:,1],depth_temps[:,2])
ts_temps=predict(predicted_temps,dr[1:50])
bottom_grad=(ts_temps[50]-ts_temps[49])/(dr[50]-dr[49])
append!(ts_temps,ts_temps[50]+bottom_grad*[20,40,60,80])



depth_age=float(cs[3:end,1:2])
ages=InterpIrregular(depth_age[:,1],depth_age[:,2],1,InterpNearest)
ts_ages=map(x->ages[x],dr)



sv=Array(Float64,0)
dages=map(i->ages[i]-ages[i-1],2:length(ages))#  ages[2:end]-ages[1:end-1]
ddepths=depth_age[2:end,2]-depth_age[1:end-1,2]


ddepthdtau=ddepths[1:3406]./dages[1:3406]
ages_good=map(x->ages[x],1:3406)
dages_good=ages_good[2:end]-ages_good[1:end-1]
#do local regression
lm=loess(ages_good,ddepthdtau)
dageslm=loess(depth_age[1:3405,2],dages_good)
smoothedDages=predict(dageslm,depth_age[1:3405,2])



pred=predict(lm,ages_good)
pred=pred/pred[1]
predictedThinning=predict(lm,ts_ages)
predictedThinning=predictedThinning/predictedThinning[1]
vertStrain=(pred[2:end]-pred[1:end-1])./smoothedDages
svslm=loess(ages_good[1:end-1],vertStrain)
smoothedVertStrain=(pred[1:end-1]./pred[2:end])./(pred[1:end-1].*smoothedDages[1:end])
svsi=InterpIrregular(ages_good,smoothedVertStrain,1,InterpNearest)

ts_smoothedVertStrain=map(x->svsi[x],ts_ages)

function dj(x,kink_depth)
    if x< kink_depth 
        return 0
    else 
        return (x-kink_depth)/1000
    end
end
     
dt=5e-5
fab.epsdot=zeros(3,3,1)
fab.vort=zeros(3,3,1)
sv=Array(Float64,0)

#initialize
w=1.
let(i=9)
    pars.dt=dt*(ts_ages[i+1]-ts_ages[i])*10
    com=ts_smoothedVertStrain[i]
    ss=dj(dr[i],2000)*100
    fab.temp=ts_temps[i]
    fab.epsdot[1,1]=-2*com
    fab.epsdot[2,2]=com
    fab.epsdot[3,3]=com
    fab.epsdot=-fab.epsdot
end
s=zeros(3)
count=0
while (w>0.187183) & (count<100)
    count+=1
    fabE(pars,fab,jefferysRHS)
    s=svd(fab.p[:,:,1])[2]
    w=min(s)/norm(s)
    print(w,"\n")
end
pout=Array(Float64,3,length(fab.p[1,:]),54)
  for i=1:length(ts_ages)-1
    pars.dt=dt*(ts_ages[i+1]-ts_ages[i])
    com=ts_smoothedVertStrain[i]
    ss=dj(dr[i],2000)*5
    fab.temp=ts_temps[i]
    fab.epsdot[3,1]=ss
    fab.epsdot[1,3]=ss
    fab.epsdot[1,1]=-2.*com
    fab.epsdot[2,2]=com
    fab.epsdot[3,3]=com
    fab.vort[3,1]=ss
    fab.vort[1,3]=-ss
    fab.epsdot=-fab.epsdot
    fabE(pars,fab,jefferysRHS)
    pout[:,:,i]=p[:,:,1]
    append!(sv,sort(svd(fab.p[:,:,1])[2]))
    end
sv=reshape(sv,(3,int(length(sv)/3)))
schmidtPlot(fab.p);plt.show()
for i=1:size(sv)[2]
  sv[:,i]=sv[:,i]/norm(sv[:,i])
end

plt.plot(ts_ages[2:54],sv[1,:]');plt.plot(ts_ages[2:54],x[1,2:54]');plt.show()

smoothedLayerThickness=ddepthdtau

function thinning(dages,ddepths)
    return (ddepths[3:end]+ddepths[1:end-2]-2*ddepths[2:end-1])/(dages[2:end-1])
  end

#get thinning function -> strain rate

    
       
#get singular values and directions

