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
  svs[d]=x[1]/norm(x)
  cd("..")
  end
cd("..")

#now interpolate depth-age
dr=sort(float(dr))

cs=readdlm("timescale.csv",',')
cd("..")


depth_age=float(cs[3:end,1:2])
ages=InterpIrregular(depth_age[:,1],depth_age[:,2],1,InterpNearest)
ts_ages=ages[dr]


sv=Array(Float64,0)
dages=ts_ages[2:end]-ts_ages[1:end-1]
  for i=1:length(dages)
    pars.dt=dt*dages[i]
    fabE(pars,fab,jefferysRHS)
    x=(svd(fab.p[:,:,1]))[2]
    sv=append!(sv,[max(x)/norm(x)])

    end

plt.plot(ts_ages[2:end],svs[2:end])
plt.plot(ts_ages[2:end],sv)
plt.show()

function wrapper(ages,p,ts_svs,fab,pars,jefferysRHS)
  t_c=p[1]
  fab.epsdot[1,1]=p[2]
  fab.epsdot[2,2]=p[3]
  fab.epsdot[3,3]=p[4]
  fab.epsdot[2,3]=p[5]
  fab.epsdot[1,3]=p[6]
  fab.epsdot[1,2]=p[7]
  symmetrize!(fab.epsdot)
  fab.vort[3,2]=-p[8]
  fab.vort[2,3]=p[8]
  fab.epsdot[2,3]=p[8]
  fab.epsdot[3,2]=p[8]
  function objective(ages,p)
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
      
       
#get singular values and directions

