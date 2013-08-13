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

function imq(r,ep=1::Number)
  return 1/sqrt(1+r.^2*ep.^2)
  end

function lapIMQ(r,ep=1::Number)
  return ep^.2*(ep.^2*r.^2-2)./(ep.^2*r.^2+1).^2.5
  end

function d2imq(x,x0,i::Int,j::Int,eps::Number)
  r=x-x0
  r2=(r*r')[1]
  if i !=j
    return 3*eps^2*r[i]*r[j]/(eps*r2+1)^2.5
    else
      return eps*(eps*(3*r2[i]-sum(r2))-1)/(eps*r2+1)^2.5
    end
  end

function lan(x0,x,C,ep)
  a=[:1,:2,:3,:3,:3,:2]
  for k=1:6
    u[k]=:(d2imq(x,x0,$(a[$k]),j,ep))
    end
  for k=1:6
    d[k]=:(sum(C[:,$k].*u))
    end
  for i=1:3
    L
  for i=1:6
    u[i]=d2imq(x,x0,i,

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
function getWeights(kd::PyObject,coors,L::Function,bnd_index::Int,n::Int,nnn::Int)
  (d,inds0)=kd[:query](coors[1:bnd_index-1,:],nnn)
  inds=inds0+1
  w=Array(Float64,bnd_index-1,nnn)
  S=Array(Float64,nnn,nnn)
  for i=1:bnd_index-1
    #generate the weights matrix
    S=Float64[imq(d[i,:],0.1)]
    Lh=L(d[i,:])
    #generate the augmented matrix
    S=[S ones(nnn)
       ones(nnn)' 0]
    Lh=[Lh',0]
    w[i,:]=(S\Lh)[1:nnn]
    end
  return (w,d,inds)
  end

#getweights for 3 dep vars 
function getWeights(kd::PyObject,coors,C,L::Function,bnd_index::Int,n::Int,nnn::Int)
  (d,spinds0)=kd[:query](coors[1:bnd_index-1,:],nnn)
  spinds=spinds0+1
  #Move from spatial to spatial/component indexing
  #i where i mod 3 == 1 is x component of spinds[i], ' '==2 is y, etc
  bin=3*bnd_index
  nnnc=3*nnn

  inds=Array(Float64,bin-3,nnnc)
  w=Array(Float64,bin-3,nnnc)
  S=Array(Float64,nnnc,nnnc)
  for i=1:bnd_index-1
    #generate the weights smatrix
    S=Float64[imq(d[i,:],0.1)]
    Lh=L(d[i,:])
    #generate the augmented matrix
    S=[S ones(nnn)
       ones(nnn)' 0]
    Lh=[Lh',0]
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
