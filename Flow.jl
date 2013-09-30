module Flow
export FlowParams,Lfl,imq,d2imq,genrRHS,genrSystem
using PyCall
@pyimport scipy.spatial as sp

type FlowParams
  coors::Array{Float64,2}
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
#[u_z]  [:
#[p]    [
#Define radial basis function and second derivatives
function imq(x,x0,ep=1::Number)
  r2=sum([x-x0].*[x-x0])
  return (1/sqrt(1+ep^2.*r2))
  end

#function imq(r,ep=1::Number)
#  return 1/sqrt(1+r.^2*ep.^2)
#  end

function lapIMQ(r,ep=1::Number)
  return ep^.2*(ep.^2*r.^2-2)./(ep.^2*r.^2+1).^2.5
  end

function dimq(x,x0,i,ep)
  r2=sum([x-x0].*[x-x0])
  return -ep^2*(x[i]-x0[i])*(1/sqrt(1+ep^2*r2))
  end

function d2imq(x,x0,i::Int,j::Int,ep::Number)
  r=x-x0
  r2=(r*r')[1]
  if i !=j
    return 3*ep^2*r[i]*r[j]/(ep*r2+1)^2.5
    else
      return ep*(ep*(3*r2-sum(r2))-1)/(ep*r2+1)^2.5
    end
  end

function Lfl(x::AbstractArray,x0::AbstractArray,C::AbstractArray,ep::Number)
  l=zeros(3)
  for i=([1,1],[2,2],[3,3],[2,3],[1,3],[1,2])
    l[1]=sum((C[i,1]+C[i,6]+C[i,5])*d2imq(x,x0,i[1],i[2],ep))+l[1]
    l[2]=sum((C[i,6]+C[i,2]+C[i,4])*d2imq(x,x0,i[1],i[2],ep))+l[2]
    l[3]=sum((C[i,5]+C[i,4]+C[i,3])*d2imq(x,x0,i[1],i[2],ep))+l[3]
    end
  return l
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
  n=length(fpar.coors[1,1]) #Make sure the first index is site #
  kd=sp.cKDTree(fpar.coors)
  (w,J,inds)=getWeights(kd,fpar.coors,fpar.L,fpar.bnd_index,n,fpar.nnn)
  I=inds'[:]
  #J=Array(Int,length(I))
  (Ibc,Jbc,Vbc)=applyBC(fpar.bnd_index,n) #
  for i=1:fpar.bnd_index-1
    J[(i-1)*fpar.nnn+1:i*fpar.nnn]=i
    end
  J=[J,Jbc]
  I=[I,Ibc]
  V=[w,Vbc]
  return sparse(J,I,V)
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
    S=[imq(d[i,:],0.1)]
    Lh=L(d[i,:])
    #generate the augmented matrix
    S=[S ones(nnn)
       ones(nnn)' 0]
    Lh=[Lh',0]
    w[i,:]=(S\Lh)[1:nnn]
    w=reshape(w',length(w))
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
  nnnt=3*nnnc #length of weights vector for three eqns at single starred point
  inds=Array(Float64,(bin-3)*nnnt)
  w=Array(Float64,(bin-3)*nnnt)
  w2=Array(Float64,nnn,6) #The weights for each du^2/dx_i^2
  J=Array(Float64,(bin-1)*nnnt)
  Lh=Array(Float64,nnnc+1)
  Lh[end]=0 #For augmented system
#convert to array indexing from spatial indexing
#inds[(i-1)*nnnc+1:i*nnnc]=[spinds,spinds+1,spinds+2]
  for i=1:(bin-1)
    #inds[(i-1)*nnnt
    J[(i-1)*nnnt+1:(i-1)*nnnt+nnnc]=(i-1)*3+1
    J[(i-1)*nnnt+nnnc+1:(i-1)*nnnt+2*nnnc]=(i-1)*3+2
    J[(i-1)*nnnt+2*nnnc+1:(i-1)*nnnt+3*nnnc]=(i-1)*3+3
    end
#match row indices
#S=Array(Float64,nnnc,nnnc)
#S is  phi(x1-x1):phi(x1-x[1:nnn]) phi(x1-x[1:nnn]) phi(x1-x[1:nnn])
#      ...           ...        ...
#      phi(xn-x1) ...
#For weights, where first set is for u, then v, then w
  for i=1:bnd_index-1
    #get the nonzero column indices for rows 3*i-1,3*i,3*i+1
    #first eqn
    inds[(i-1)*nnnt+1:3:(i-1)*nnnt+nnnc]=3*spinds[i,:]
    inds[(i-1)*nnnt+1:3:(i-1)*nnnt+nnnc]=3*spinds[i,:]+1
    inds[(i-1)*nnnt+1:3:(i-1)*nnnt+nnnc]=3*spinds[i,:]*2
    #second
    inds[(i-1)*nnnt+nnnc+1:3:(i-1)*nnnt+2*nnnc]=3*spinds[i,:]
    inds[(i-1)*nnnt+nnnc+1:3:(i-1)*nnnt+2*nnnc]=3*spinds[i,:]+1
    inds[(i-1)*nnnt+nnnc+1:3:(i-1)*nnnt+2*nnnc]=3*spinds[i,:]*2
    #third
    inds[(i-1)*nnnt+2*nnnc+1:3:(i-1)*nnnt+3*nnnc]=3*spinds[i,:]
    inds[(i-1)*nnnt+2*nnnc+1:3:(i-1)*nnnt+3*nnnc]=3*spinds[i,:]+1
    inds[(i-1)*nnnt+2*nnnc+1:3:(i-1)*nnnt+3*nnnc]=3*spinds[i,:]*2

    #generate the weights matrix
    #Needs casting to avoid Array{Any...}
    for j=1:nnn
      for k=1:nnn
        S[j,k]=imq(coors[spinds[i,j],:],coors[spinds[i,k],:])
        end
      end
      S=[S ones(nnn)
        ones(nnn) 0]
    let d2=Arrt(Float64,nnn+1);dv=[[1,1],[2,2],[3,3],[2,3],[1,3],[1,2];]
        S1=Array(Float64,nnn+1,nnn+1)
      d2[7]=0
      for j=1:nnn
        for k=1:nnn
          S1[j,k]=imq(coors[spinds[i,j],:],coors[spinds[i,k],:])
          end
        end
      S1=[S1 ones(nnn)
      ones(nnn) 0]
      for j=1:6
        for k=1:nnn
          d2[k]=d2imq(coors[spinds[i,k],:],coors[spinds[i,1],:],dv[j,1],dv[j,2],ep)      
          end
          w2[:,j]=(S1\d2)[1:nnn]
        end
      end #let
    #S1=Float64[imq(coors[j,:]-coors[k,:]) for j in spinds[i,:],k in spinds[i,:]] 
    #w[i,j] is the weights to approx the second derivative i=1:6 (voigt) at spatial
    #point j. Now, w[(i-1)*nnnt+j:3:i*nnnt] is weights for j=u,v,w velocity components for <<i'th spatial point>> and equation << 1 >> (of 3 coupled)
    let uc=C[6,:]+C[5,:]+C[1,:]
      #first equation S[x,x],x+S[x,y],x+S[x,z],x=whatever
      w[(i-1)*nnnt+1:3:(i-1)*nnnt+nnnc]=uc[1]*wc[1]+0.5*uc[5]*wc[:,5]+0.5*uc[6]*wc[:,6]
      w[(i-1)*nnnt+2:3:i*nnnc]=uc[2]*wc[:,6]+0.5*uc[4]*wc[:,6]+0.5*uc[6]*wc[:,1]
      w[(i-1)*nnnt+1:3:i*nnnc]=uc[3]*wc[:,5]+0.5*uc[:,4]*wc[:,6]+0.5*uc[5]*wc[:,1]
      #second eqn
      uc=C[6,:]+C[2,:]+C[5,:]
      w[(i-1)*nnnt+nnnc+1:3:(i-1)*nnnt+2*nnnc]=uc[1]*wc[:,6]+0.5*uc[5]*wc[:,6]+0.5*uc[6]*wc[:,1]
      w[(i-1)*nnnt+nnnc+2:3:(i-1)*nnnt+2*nnnc]=uc[2]*wc[:,2]+0.5*uc[4]*wc[:,5]+0.5*uc[6]*wc[:,6]
      w[(i-1)*nnnt+nnnc+3:3:(i-1)*nnnt+2*nnnc]=uc[3]*wc[:,4]+0.5*uc[4]*wc[:,3]+0.5*uc[5]*wc[:,6]
      #third one
      uc=C[5,:]+C[4,:]+C[3,:]
      w[(i-1)*nnnt+2*nnnc+1:3:i*nnnt]=uc[1]*wc[:,5]+0.5*uc[5]*wc[:,3]+0.5*uc[6]*wc[:,4]
      w[(i-1)*nnnt+2*nnnc+2:3:i*nnnt]=uc[2]*wc[:,4]+0.5*uc[4]*wc[:,3]+0.5*uc[6]*wc[:,5]
      w[(i-1)*nnnt+2*nnnc+3:3:i*nnnt]=uc[3]*wc[:,3]+0.5*uc[4]*wc[:,4]+0.5*uc[5]*wc[:,5]
      end

#function Lfl(x::AbstractArray,x0::AbstractArray,C::AbstractArray,ep::Number)
#  l=zeros(3)
#  for i=([1,1],[2,2],[3,3],[2,3],[1,3],[1,2])
#    l[1]=sum((C[i,1]+C[i,6]+C[i,5])*d2imq(x,x0,i[1],i[2],ep))+l[1]
#    l[2]=sum((C[i,6]+C[i,2]+C[i,4])*d2imq(x,x0,i[1],i[2],ep))+l[2]
#   1 l[3]=sum((C[i,5]+C[i,4]+C[i,3])*d2imq(x,x0,i[1],i[2],ep))+l[3]
#    end
#  return l
#  end
    #vc[:]=sum((C[i,6]+C[i,2]+C[i,4])*d2imq(x,x0,i[1],i[2],ep))+l[2]
    #l[3]=sum((C[i,5]+C[i,4]+C[i,3])*d2imq(x,x0,i[1],i[2],ep))+l[3]
    end
    #don't fuck up the indexing when calling.
    return (w,J,inds)
  end #function

#for dirichletBCs
function applyBC(bnd_index,n)
  I=[bnd_index:n]
  J=I
  V=ones(n-bnd_index+1)
  return (I,J,V)
  end
end #module
