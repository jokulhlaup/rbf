module Plotting
using jefferys,PyCall
@pyimport matplotlib.pyplot as plt
export rextalPlot,schmidtPlot

#converts from normals to a projection in polar coordinates.
function normal2polar(c)
#Expects the first index to be length three.
  n=length(c)
  m=convert(Int,length(c)/3)
  rem(n,3)==0?nothing:error("Length of c array ain't divisible by three.")
  theta=Array(Float64,convert(Int,n/3))
  r=Array(Float64,convert(Int,n/3))
  for i=1:m
    r[i]=sqrt(c[3*i-2]^2+c[3*i-1]^2)
    theta[i]=atan2(c[3*i-2],c[3*i-1])
    end
  return (r,theta)
  end

function schmidtPlot(c)
  (r,theta)=normal2polar(c)
  ax=plt.subplot(111, polar=true)
  p=plt.scatter(theta,r)
  return p
  end


function eqArea(c)
    n = size(c,2)
    for i=1:n
        if c[3,i] <0
            c[:,i]=-c[:,i]
        end
    end
    (X,Y) = (sqrt((1+c[3,:])/2).*2.*c[1,:], sqrt((1+c[3,:])/2).*2.*c[2,:])
    r=sqrt(X.*X .+ Y.*Y)

    theta = atan2(Y,X)
    ax=plt.subplot(111,polar=true)
    p=plt.scatter(theta,r)
    return p
end

#function kdE(
function rextalPlot(c,r)
  n=length(r)
  cn=filter(i->r[i]>1e-6,1:n)  
  return schmidtPlot(c)
  end
end
