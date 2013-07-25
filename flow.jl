module Flow
export imq,d2imq
using PyCall
@pyimport scipy.spatial as sp


#Need to get a 6x6 viscosity matrix C
#routine to build sparse matrix by COLUMN 
#(julia uses CSC storage)
#to solve 
#[u_x]' [
#[u_y]  [
#[u_z]  [
#[p]    [
#Define radial basis function and second derivatives
function imq(x::Array{Float64,1},x0::Array{Float64,1})
  r2=(x-x0)'*(x-x0)
  return (1/sqrt(1+eps^2*r2))
  end
function d2imq(x,x0,eps,i,j)
  r=x-x0
  r2=r'r
  if i !=j
    return 3*eps^2*r[i]*r[j]/(eps*r2+1)^2.5
    else
      return eps*(eps*(3*r2[i]-sum(r2))-1)/(eps*r2+1)^2.5
    end
  end

end #module
