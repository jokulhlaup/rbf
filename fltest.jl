 
 ##################3
#######Tests
coors=unifmesh([1:5],[1:5])
nnn=5
bnd_index=25
L=(x0,x1)->d2imq(x1,x0,1,1,0.1)+d2imq(x1,x0,1,1,0.1) 
fpar=FlowParams(coors,nnn,bnd_index,L)   

A=genrSystem(fpar)
