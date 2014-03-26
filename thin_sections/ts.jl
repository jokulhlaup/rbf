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
p=Dict()
svs=Array(Float64,length(dr))
for d=1:length(dr)
  cd(dr[d])
  rf=readdlm("c-axes.txt",'\n')
  i=int(dr[d])
  c[i]=readif(rf)
  #convert to rads
  c[i]=c[i]*pi/180
  p[i]=Array(Float64,3,length(c[i][1,:]))
  p[i][1,:]=cos(c[i][1,:])
  p[i][2,:]=sin(c[i][1,:]).*cos(c[i][2,:])
  p[i][3,:]=sin(c[i][1,:]).*sin(c[i][2,:])
  x=svd(p[i])[2]
  svs[d]=minimum(x)/norm(x)
  cd("..")
  end

cd("..")
cs=readdlm("timescale.csv",',')

#now interpolate depth-age
dr=sort(float(dr))



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

fab.epsdot=zeros(3,3,1)
  for i=1:length(ts_ages)-1
    com=ts_smoothedVertStrain[i]
    fab.epsdot[3,3]=-com
    fab.epsdot[2,2]=-com
    fab.epsdot[1,1]=2*com
    pars.dt=dt*(ts_ages[i+1]-ts_ages[i])
    fabE(pars,fab,jefferysRHS)
    append!(sv,[minimum(svd(fab.p[:,:,1])[2])])
    end



smoothedLayerThickness=ddepthdtau

function thinning(dages,ddepths)
    return (ddepths[3:end]+ddepths[1:end-2]-2*ddepths[2:end-1])/(dages[2:end-1])
  end

#get thinning function -> strain rate

    
       
#get singular values and directions

