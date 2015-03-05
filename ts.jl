require("lsap/Assignment.jl")
require("Constructors.jl")
using Loess,Grid
@pyimport matplotlib.pyplot as plt

function emdSM(p)
    Utils.proj2UpHem!(p)
    n=size(p,2)
    a=p'*[0,0,1]
    a[abs(a) .> 1.]=0.99999
    return sum(abs(acos(a)))/n 
    end



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

cd("thin_sections/C-axisdatabase2/sections/")
dr=readdir()
c=Dict()
pd=Dict()
svs=Array(Float64,length(dr))
emdss_core=Array(Float64,length(dr))
x=Array(Float64,3,length(dr))
for d=1:length(dr)
  cd(dr[d])
  istr=dr[d]
  rf=readdlm("$istr c-axes.txt",'\n')
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
    pd[i][3,j] = -pd[i][3,j] ;
      end
    end
  x[:,d]=eigvals(pd[i]*pd[i]')/size(pd[i],2)
  svs[d]=minimum(x[:,d])/norm(x[:,d])
  emdss_core[d]=emdSM(pd[i])
  cd("..")
  end

cd("../..")
cs=readdlm("timescale.csv",',')

for d=1:length(dr)
    i=int(dr[d])
    x[:,d]=eigvals(pd[i]*pd[i]')/size(pd[i],2)
end


#now interpolate depth-age
dr=sort(float(dr))
ns=length(dr)
depth_temps=readdlm("temp.TXT")
predicted_temps=loess(depth_temps[:,1],depth_temps[:,2])
ts_temps=predict(predicted_temps,dr[1:end-4])
bottom_grad=(ts_temps[end]-ts_temps[end-1])/(dr[79]-dr[78])
append!(ts_temps,ts_temps[end]+bottom_grad*[20,40,60,80])



depth_age=convert(Array{Float64,2},cs[3:end,1:2])
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
        return (x-kink_depth)/300
    end
end
     
(fabE,fab,pars)=Constructors.mkFab(300)
dt=2e-5
fab.epsdot=zeros(3,3,1)
fab.vort=zeros(3,3,1)
function init2svd(fab,smsvd)
  #initialize
  w=1.
  let(i=1)
    pars.dt=dt*(ts_ages[i+1]-ts_ages[i])*10
    com=ts_smoothedVertStrain[i]
    ss=dj(dr[i],2000)*1000
    fab.temp=ts_temps[i]+10
    fab.epsdot[1,1]=2.0*com
    fab.epsdot[2,2]=0.0
    fab.epsdot[3,3]=-2.0*com
    end
  s=zeros(3)
  count=0
  while (w>smsvd) & (count<100)
    count+=1
    fabE(pars,fab,jefferysRHS)
    s=svd(fab.p[:,:,1])[2]
    w=minimum(s)/norm(s)
    print(w,"\n")
    end
  end

function getGirdle(n)
  thetas=linspace(-pi/2,pi/2,n)
  return [zeros(n)',sin(thetas)',cos(thetas)']
  end

jefferysRHS=jefferys.jefferysRHS
init2svd(fab,0.8)

function emdSM(p)
    Utils.proj2UpHem!(p)
    n=size(p,2)
    a=p'*[0,0,1]
    a[abs(a) .> 1.]=0.99999
    return sum(abs(acos(a)))/n 
    end

for i=1:length(dr)
    n=size(pd[int(dr[i])],2)
    girdle=getGirdle(n)
    pd[int(dr[i])]=Utils.alignFabrics(girdle,pd[int(dr[i])])[3]

end

function evThruCore(p,get_r)
emdss=Array(Float64,ns)
emdg=Array(Float64,ns)
svm=Array(Float64,3,ns)
n=size(p)[2]

rs=Array(Float64,n,ns)
dt=2e-5
sma=repmat([0. 0. 1.],n)'
(fabE,fab,pars)=Constructors.mkFab(n)
fab.r=map(x->rand()<0.5?x=0:x,fab.r)
fab.p[:,:,1]=p;
emdist=zeros(ns)
grmobs=ones(ns)
#grmobs[21:end]=3.
girdle=getGirdle(n)
pout=zeros(3,n,ns)
fab.vort=zeros(3,3,1)
pars.nrk=100
  for i=1:ns-1
    print("i = $i !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!")
    pars.dt=dt*(ts_ages[i+1]-ts_ages[i])
    com=2*ts_smoothedVertStrain[i]
    ss=-dj(dr[i],2383)*2
    fab.temp=ts_temps[i]+ 10
    println(fab.temp, "temps")
    fab.epsdot[3,1]=ss
    fab.epsdot[1,3]=ss
    fab.epsdot[1,1]=com
    fab.epsdot[2,2]=-0.3*com
    fab.epsdot[3,3]=-0.7*com
    fab.vort[1,3]=ss
    print(fab.epsdot)
    fab.vort[3,1]=-ss
#    fab.vort=zeros(size(fab.vort))
    #fab.vort=fa-fab.epsdotb.vort ##
    #fab.epsdot=-fab.epsdot ##
    #fab.epsdot[1,3]= -fab.epsdot[1,3];fab.epsdot[3,1]= -fab.epsdot[3,1]
    fab.grmob=grmobs[i]
#    println("1 ",fab.epsdot)
    fabE(pars,fab,jefferysRHS)
#    println("2 ",fab.epsdot)
    pout[:,:,i]=fab.p[:,:,1]
    println(fab.str[10,1]) 
    println(pars.dt, " dt")
    emdss[i]=emdSM(fab.p[:,:,1])
    rs[:,i]=fab.r[:,1]
    global rad=fab.r
    global strain=fab.str
    emdg[i]=Utils.earthMoversDist(fab.p[:,:,1],girdle)[1]/n
    foo=eigvals(fab.p[:,:,1]*fab.p[:,:,1]')
    print(typeof(foo))
    svm[:,i]=sort(foo)
    #svm[i]=minimum(foo)/norm(foo)
    println(emdss[i])
#    append!(sv,sort(svd(fab.p[:,:,1])[2]))
    println("\n step",i)
    end
  #sv=reshape(sv,(3,int(length(sv)/3)))  
#  for i=1:size(ddsv)[2]
#    sv[:,i]=sv[:,i]/norm(sv[:,i])
#    end
 println(typeof(pout),typeof(emdss),typeof(emdg),typeof(fab))
  return (pout,emdss,emdg,fab,rs,svm)
end



for d=1:length(dr)
    i=int(dr[d])
    x[:,d]=eigvals(pd[i]*pd[i]')/size(pd[i],2)
end


function get_mean_rs(comp_out)
    rs=comp_out[end-1]
    (ng,nc,ns)=size(rs)
    rmean=zeros(nc,ns)
    rmg=zeros(nc)
    for j=1:nc
        for i=1:ns
            rmean[j,i]=mean(filter(x->x>1e-6,rs[:,j,i]))
        end
        rmg[j]=mean(rs[j,:])
    end
    
    return rmean,rmg
end

function nzeigs(comp_out,nrun=1)
    eignz=zeros(3,83);
    for i=1:83
       bark=Utils.filterZeros(comp_out[1][:,:,i,nrun],comp_out[end-1][:,nrun],1);
       nnz=size(bark,2)
       eignz[:,i]=eigvals(bark*bark')/nnz;
   end
   return eignz
end


function plot_eigs(comp_out)
    nrun=size(comp_out[1])[end]
    for i=1:nrun
        eignz=nzeigs(comp_out,i)
        plt.plot(eignz[3,:]')
    end
    plt.plot(x[3,:]',linestyle="-.")
    plt.show()
end

evThruCore(p)=evThruCore(p,true)[1:4]
function trym1(m,n,get_r)
    svm=Array(Float64,3,ns,m) 
    emdg=Array(Float64,ns,m)
    emdss=Array(Float64,ns,m)
    pout=zeros(3,n,ns,m)
    r=Array(Float64,n,ns,m)
    i=1
    while (i <= m)
       inds=rand(1:n,n)
       p=pd[int(dr[1])][:,inds]
       inds=rand(1:size(p,2),n)
       try 
           (pout[:,:,:,i],emdg[:,i],emdss[:,i],fab,r[:,:,i],svm[:,:,i])=evThruCore(p,true)
           i+=1
       catch
           print("THERE WAS AN ERROR")
           sleep(5)
           continue
       end
    end
   return (pout,emdg,emdss,fab,r,svm)
   end

trym1(m,n)=trym1(m,n,true)[1:4]
function plotm(m,emdg)
    for i=1:m-1
        plt.plot(emdg[1:53,i])
    end
    plt.show()
end

#function getsvs(pout)
#    svss=Array(Float64,size(pout,3),size(pout,4))
#    svg=deepcopy(svss)
#    for i=1:size(pout,4)
#        for j=1:size(pout,3)
#            svss[i,j]=


function getquantile(mat, p)
    m,n=size(mat)
    q=zeros(m)
    for i=1:m
        q[i]=quantile(vec(mat[i,:]),p)
    end
    return q
end
function plot_eig_quants(comp_out,eig_num=3)
    (th,n,nc,ns)=size(comp_out[1])
    eig_mat=zeros(nc,ns)
    for i=1:ns
        eig_mat[:,i]=nzeigs(comp_out,i)[3,:]
    end
    pl=pyQuants(eig_mat[1:nc-1,:],[0.025,0.5,0.975],["asdf","median"])
    return pl
end

function pyQuants(mat,ps,labs)
    fig,ax=plt.subplots(1)
    q=Array(Float64,82,length(ps))
    for i=1:3
       q[:,i]=getquantile(mat,ps[i])
    end
    ax[:plot](ts_ages[1:82]/1000,q[:,2],label=labs[2])
    ax[:fill_between](ts_ages[1:82]/1000,q[:,1],q[:,3],facecolor="yellow",alpha=0.5,label="95% empirical quantile range")
    ax[:plot](ts_ages[1:82]/1000,mat[:,1],label="sample a")
    ax[:plot](ts_ages[1:82]/1000,mat[:,2],label="sample b")
    ax[:plot](ts_ages[1:82]/1000,x[3,1:82]',label="core")
    ax[:legend](loc="upper left")
    ax[:set_xlabel]("ice age (ka)")
    ax[:set_ylabel]("Largest eigenvalue")
        #x=ts_ages[1:53]/1000,y=q[1:53,i],label=labs[i])
        #pl_df=vcat(pl_df,df)
#        plt.plot(ts_ages[1:53]/1000,q[1:53], color="k",label=labs[i],marker=m[i]
#        append!(delayed_plot.args,[layer(x=ts_ages[1:53]/1000,y=q[1:53],Geom.line)])
   return (fig,ax)
    end
 
function plotQuants(mat,ps,labs)
    m=["+","2",","]
    delayed_plot=:(plot())
    q=Array(Float64,53,length(ps))
    pl_df=DataFrames.DataFrame()
    for i=1:length(ps)
        q[:,i]=getquantile(mat,ps[i])
        df=DataFrames.DataFrame(x=ts_ages[1:53]/1000,y=q[1:53,i],label=labs[i])
        pl_df=vcat(pl_df,df)
#        plt.plot(ts_ages[1:53]/1000,q[1:53], color="k",label=labs[i],marker=m[i]
#        append!(delayed_plot.args,[layer(x=ts_ages[1:53]/1000,y=q[1:53],Geom.line)])
    end
    pl_df=vcat(pl_df,DataFrames.DataFrame(x=ts_ages[1:53]/1000,y=mat[1:53,1],label="sample 1"))
  #  pl_df=vcat(pl_df,DataFrame(x=ts_ages[1:53]/1000,y=mat[1:53,2],label="sample 2"))
  #  pl_df=vcat(pl_df,DataFrame(x=ts_ages[1:53]/1000,y=emdss_core[1:53],label="WAIS"))
#     append!(delayed_plot.args,[layer(x=ts_ages[1:53]/1000,y=emdg[1:53,1],Geom.line)])
#     append!(delayed_plot.args,[layer(x=ts_ages[1:53]/1000,y=emdg[1:53,2],Geom.line)])
#    eval(delayed_plot)
#    plt.plot(ts_ages[1:53]/1000,emdg[1:53,10],"--",color="k",label="sample a");plt.plot(ts_ages[1:53]/1000,emdg[1:53,40],"-.",color="k", label="sample b")
#    plt.xlabel("ice age (ka)")
#    plt.ylabel("Fabric earth mover distance from girdle (dimensionless)")
#    plt.legend()
   
   plot(pl_df, x="x", y="y",color="label", Geom.line, Guide.xlabel("ice age (ka)"), Guide.ylabel("earth mover distance from girdle (dimensionless)",orientation=:vertical), Scale.discrete_color())
end

function r2(mat)
    m,n=size(mat)
    TSS=zeros(m);SSE=zeros(m)
    s=sum(a,2)/n
    mean_tot=mean(mat)
    mat2=deepcopy(mat)
    for i=1:m
        SSE[i]=sum((mat2[:,i].-s).^2)
        TSS[i]=sum((mat2[:,i]-mean(mat2[:,i])).^2)
    end
    return 1 - SSE./TSS
end
emdstss=zeros(ns)
emdsg=zeros(ns)
for i=1:40
    global pda=Utils.proj2UpHem!(pd[int(dr[i])])
    pda=Utils.rotp(pi/2,Utils.alignFabrics(girdle,pda)[3])
    n1=size(pda,2)
    girdle=getGirdle(n1)
    emdstss[i]=emdSM(pda)
    emdsg[i]=Utils.earthMoversDist(pda,girdle)[1]/n1
   # emdgtss[i]=emdg(
   # emdssp[i]=emdSM(p)
    print(i)
end
#function evDiff(m)
#    TmpAdjust=
#    sv=Array(Float64,0)
#    n=size(p)[2]
#    (fabE,fab,pars)=Constructors.mkFab(n)
#    fab.p[:,:,1]=p;
#    grmobs=ones(ns)
#    grmobs[21:end]=100
#    pout=zeros(3,n,ns)
#    fab.vort=zeros(3,3,1)


schmidtPlot(fab.p);plt.show()
for i=1:size(sv)[2]
  sv[:,i]=sv[:,i]/norm(sv[:,i])
end
function evAll(fab)
  for i=1:53
    n=size(pd[int(dr[i])],2)
    girdlequiv=getGirdle(n)
    sv[i,:,:]=evThruCore(pd[int(dr[i])])
    println(i)
    end
  return sv
  end



function sPn(n)
  schmidtPlot(pout[:,:,n])
  plt.show()
  schmidtPlot(pd[int(dr[n])])
  plt.show()
  end

plt.plot(ts_ages[2:ns],sv[1,:]');plt.plot(ts_ages[2:ns],x[1,2:ns]');plt.show()

smoothedLayerThickness=ddepthdtau

function thinning(dages,ddepths)
    return (ddepths[3:end]+ddepths[1:end-2]-2*ddepths[2:end-1])/(dages[2:end-1])
  end

#get thinning function -> strain rate

    
       
#get singular values and directions

