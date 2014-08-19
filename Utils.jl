module Utils
using Base.rand,PyCall#,Optim
include("lsap/Assignment.jl")
export voigt2Tensor,tensor2Voigt,rk4!
export halton,vdc,unifmesh,randir
export diffrandi,secondInv,binBoolInd,earthMoversDist
unshift!(PyVector(pyimport("sys")["path"]), "")
@pyimport pyutils

#function get_perms(v,depth=1)
#    n=length(v)
#    u=deepcopy(v)
#    assert(n<7, "too long")
#    perms=Array(Any,factorial(v))
#
#    if n==1
#        perms=[v]; return;
#    for i=depth+1:n #swap first and i
#        u[depth]=v[i]
#        u[i]=v[depth]
#        get_perms(u,depth+1)
        



function proj2UpHem!(p)
  n=size(p)[2]
  for i=1:n
    if p[3,i]<0.
        p[:,i]= -p[:,i]
      end
    end
  return p
  end


function getRandc(n)
  p=rand(3,n)
  for i=1:n
    p[:,i]=p[:,i]/norm(p[:,i])
    end
    return p
  end

function distMat(ps1,ps2)
    m=size(ps1,2)
    n=size(ps2,2)
    costs=Array(Float64,m,n)
    for i=1:m
        for j=1:n
            dp=abs(dot(ps1[:,i],ps2[:,j]))
            if dp <=1.
                costs[i,j]=acos(dp)
            elseif 1.<dp
                costs[i,j]=0.
            end
        end
    end
   return costs
end
   


#find the horizontal axis rotation with the optimal
#alignment between p1 and p2.
#Returns (theta,p2_rotated,distance)
function alignFabrics(p1,p2)
  svp1=svd(p1)[1][:,3]
  obj(theta)=alignFabricsObj(theta,svp1,p2)
  res=Optim.optimize(obj,[0.],method=:simulated_annealing)
  return(-res.minimum,res.f_minimum,rotp(-res.minimum[1],p2))
  end
function alignFabricsObj(thetab,svp1,p2)
  theta=thetab[1];
  p2_rot=rotp(theta,p2)
  tc=svd(p2_rot)[1][:,3]
  return acos(abs(dot(tc,svp1)))
  end

function rotp(theta,p)
  rotM=zeros(3,3);rotM[1,1]=cos(theta);rotM[2,2]=cos(theta);
  rotM[3,3]=1;rotM[1,2]=-sin(theta);rotM[2,1]=-rotM[1,2];
  p_rot=rotM*p;
  return p_rot
  end
#This uses the Munkres algorithm to find 
#the emd between ps1 and ps2, where psj[spatial coors,point index]
#Munkres assigns n1 workers (ps2) to n2 jobs (ps1)
#according to cost matrix.
#If size(ps1) != size(ps2) then one needs must deal with weights. Spse p1 bigger than p2. Then have
#wts wts1 and wts2, st sum(wts1)=sum(wts2);
#Th

function earthMoversDist(ps1::Array{Float64,2},ps2::Array{Float64,2})
  (dim,n)=size(ps1)
  costs=Array(Float64,n,n)
  for i=1:n
    for j=1:n
      dp=abs(dot(ps1[:,i],ps2[:,j]))
      if dp <=1.
        costs[i,j]=acos(dp)
        elseif 1.<dp
          costs[i,j]=0.
        end
      end
  end
  tc=assignment(costs)
  return tc
  end

function emd2sm(p)
  n=size(p)[2]
  sm=repmat([0.,0.,1.]',n)'
  return earthMoversDist(p,sm)
  end
    
function assignDiffN(srtwts1,srtwts2)
  nwts1=Array(Float64,m) #New wts1 array. Some
  # of the mass from the largest m-n moved to the on
  #ones at the end m-n elements, corresponding
  #to extra pts 1:m-n
  n=length(srtwts1);m=length(srtwts2)
  sw1=sum(srtwts1);sw2=sum(srtwts2)
  assert(n<m); #sum of wts must also equal
  dn=m-n; 
  wts=[srtwts1,srtwts1[1:m-n+1]]/sw1*sw2
  end


  
  
function filterZeros(x::Array{Float64,2},y::Array{Float64,1},ind::Int64)
  m,n=size(x)
  length(y)==n?nothing:error("Dimension mismatch")
  res=Array(Float64,0)
  ind==2?x=x':nothing
  for i=1:n
     y[i]==0?nothing:append!(res,x[:,i])
  end
  ind==2?x=x':nothing
  return reshape(res,(m,int(length(res)/m)))
  end
function rk4!(f::Function,h::Float64,n::Int64,x,vort,epsdot,m)
   for i=1:n
      k1=f(x,vort,epsdot)
      k2=f(x+k1*h/2,vort,epsdot)
      k3=f(x+k2*h/2,vort,epsdot)
      k4=f(x+k3*h/2,vort,epsdot)
      x+=(1/6)*h*(k1+2*k2+2*k3+k4)
      end
   return x
   end

function voigt2Tensor(v)
  x=zeros(3,3)
  x[1,1]=v[1];x[1,2]=v[6];x[1,3]=v[5]
  x[2,3]=v[4];x[2,2]=v[2];x[3,3]=v[3]
  return x+x'
  end
function tensor2Voigt(v)
  x=zeros(6)
  x[1]=v[1,1];x[2]=v[2,2];x[3]=v[3,3]
  x[4]=v[2,3];x[5]=v[1,3];x[6]=v[1,2]
  return x
  end
function binBoolInd(x,fn,n)
  ind=1
  for i=2:n
    if fn(x[ind],x[i])==true
      ind=i
      end
    end
  return ind
  end

function vdc(n,base)
  x,denom=0.,1.
  while n>0
    x+=n%base/(denom*=base)
    n=floor(n/base)
    end
    return x   
  end
  
function halton(n::Int,dim::Int,skip::Int=1000)
  bases=[2,3,5,7,11,13,17,19,23,29]
  xs=zeros(dim,n+skip)
  for i=1:dim
    b=bases[i]
    for j=1:n+skip
      xs[(j-1)*dim+i]=vdc(j,b)
      end
    end
  return xs[:,skip:end]
  end

function unifmesh(x,y)
  m=length(x)
  n=length(y)
  r=Array(Float64,m*n,2)
  for i=1:n
    r[m*(i-1)+1:i*m,:]=[x fill(y[i],m)]
    end
    return r
  end

#Get the second invariant of a 3x3 symmetric tensor in Voigt notation.
function secondInv(x::Array{Float64,1})
  return x[1]*x[2]+x[1]*x[3]+x[2]*x[3]-x[4]^2-x[5]^2-x[6]^2
  end

secondInv(x::Array{Float64,2})=secondInv([x[1,1],x[2,2],x[3,3],x[2,3],x[1,3],x[1,2]])

function randir(dims,lo::Int64,hi::Int64)
  x=rand(Int64,dims)
  x=mod(x,hi-lo+1)+lo
  end

#creates an array of random Ints in range [lo,hi] inclusive. 
function randi(len::Int,lo::Int,hi::Int)
  x=mod(rand(Int64,len),hi-lo+1)+lo
  end
#scalar version
function randi(lo::Int,hi::Int)
  x=mod(rand(Int64),hi-lo+1)+lo
  end
#choose an Int in the range lo high but not this.
function diffrandi(this,lo,hi)
  x=this
  while x==this
    x=randi(lo,hi)
    end
  return x
  end
function randir(lo::Int64,hi::Int64)
  x=rand(Int64)
  x=mod(x,hi-lo+1)+lo
  end

function diffrandi(this,lo,hi)
  x=this
  while x==this
    x=randir(lo,hi)
    end
  return x
  end
  
 


end #module
