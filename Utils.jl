module Utils
using Base.rand
export halton,vdc,unifmesh,randir,diffrandi,secondInv

function vdc(n,base)
  x,denom=0,1
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
