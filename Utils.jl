module Utils
export halton,vdc,unifmesh

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

function secondInv(x::Array{Number,1})
  return x[1]*x[2]+x[1]*x[3]+x[2]*x[3]-x[4]^2-x[5]^2-x[6]^2

end #module
