using Utils
using Flow
using PyCall
@pyimport scipy.spatial as sp
 ##################3
#######Tests
coors=unifmesh([1:5],[1:5])
coors=rand(40,3)
nnn=5
bnd_index=25
C=rand(6,6,10)
L=(x0,x1)->d2imq(x1,x0,1,1,0.1)+d2imq(x1,x0,1,1,0.1) 
fpar=FlowParams(coors,nnn,bnd_index,Lfl)   

A=genrSystem(fpar)
