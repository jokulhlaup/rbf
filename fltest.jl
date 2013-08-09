 
 ##################3
#######Tests
using Flow
using Utils
using PyCall,Flow
@pyimport scipy.spatial as sp
coors=unifmesh([0.01:0.01:1],[0.01:0.01:1])

nnn=5
bnd_index=length(coors[:,1])
bcnodes=
L=(x0,x1)->d2imq(x1,x0,1,1,0.1)+d2imq(x1,x0,2,2,0.1) 
fpar=Flow.FlowParams(coors,nnn,bnd_index,L)   

A=Flow.genrSystem(fpar)

function d2imq(x,x0,i,j,ep=1)
  x=x-x0
  r=x'*x
  if i=j
    return 3*ep^4*x[i]^2/(ep^2*r+1)^2.5-ep^2/(ep^2*r+1)^1.5
    else
      return 0
    end
  end
    


