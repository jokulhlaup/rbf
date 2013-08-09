module Flow
export imq,d2imq,genrRHS,genrSystem
using PyCall
@pyimport scipy.spatial as sp

type FlowParams
  coors::Array{Number,2}
  nnn::Union(Array{Int,1},Int)
  bnd_index::Int
  L::Function
  end
    
#Need a dict of FabricPt.
#extract 

#Need to get a 6x6 viscosity matrix C
#routine to build sparse matrix by COLUMN 
#(julia uses CSC storage)
#to solve 
#[u_x]' [
#[u_y]  [
#[u_z]  [
#[p]    [
#Define radial basis function and second derivatives
function imq(x,x0,ep=1::Number)
  r2=sum([x-x0].*[x-x0])
  return (1/sqrt(1+ep^2*r2))
  end

function d2imq(x,x0,i::Int,j::Int,ep=1)
  x=x-x0
  r=sum(x.*x)
  de=ep.^2.*r+1
  if i==j
    return (3*ep^4*x[i].^2./de.^2.5-ep.^2./de.^1.5)[1]
    else
      return (3.*ep.^4.*x[i].*x[j]./de.^2.5)[1]
    end
  end
 

function creatDict(xs) #input of list of points
  n=length(xs[1,:])
  Xd=Dict{int,FabricPt}
  for i=1:n
    Xd[i]=FabricPt(coors)
    end
  end


function genrRHS(coors,fpar::FlowParams)
  return sin(coors[:,1])*cos(coors[:,2])
  end
#This gets the CSC matrix for the transposed system
#x'A'=b'
function genrSystem(fpar::FlowParams)
  n=length(fpar.coors[:,1]) #Make sure the first index is site #
  kd=sp.cKDTree(fpar.coors)
  (w,d,inds)=getWeights(kd,fpar.coors,fpar.L,fpar.bnd_index,n,fpar.nnn)
  w=reshape(w,length(w))
  I=inds'[:]
  J=Array(Int,length(I))
  (Ibc,Jbc,Vbc)=applyBC(fpar.bnd_index,n) #
  for i=1:fpar.bnd_index-1
    J[(i-1)*fpar.nnn+1:i*fpar.nnn]=i
    end

  J=[J,Jbc]
  I=[I,Ibc]
  V=[w,Vbc]
  return sparse(I,J,V)
  end

#d is a n-length 1d array of distances,
#where a[i]=dist(starred pt, pt i)
function getWeights(kd,coors,L::Function,bnd_index::Int,n::Int,nnn::Int)
  (d,inds0)=kd[:query](coors[1:bnd_index-1,:],nnn)
  inds=inds0+1
  w=Array(Float64,bnd_index-1,nnn)
  S=Array(Float64,nnn,nnn)
  for i=1:bnd_index-1
    #generate the weights matrix
    S=Float64[imq(coors[j,:],coors[k,:],0.1) for j in inds[i,:],k in inds[i,:]]
    Lh=Float64[L(coors[j,:],coors[i,:]) for j in inds[i,:]]
    #generate the augmented matrix
    S=[S ones(nnn)
       ones(nnn)' 0]
    Lh=[Lh,0]
    w[i,:]=(S\Lh)[1:nnn]
    end
  return (w,d,inds)
  end


#for dirichletBCs
function applyBC(bnd_index,n)
  I=[bnd_index:n]
  J=I
  V=ones(n-bnd_index+1)
  return (I,J,V)
  end
end #module
