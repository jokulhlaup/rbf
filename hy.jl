using Utils,Plotting

p0=[0,0,1]
p1=rand(3)
p0=p0/norm(p0);p1=p1/norm(p1)
p=[p0 p1]
b=0.1
sigma=float([-1/2,-1/2,1,0,0,0])
epsdot=getepsdot(p,sigma,b)
  
function evFabric(p,epsdot,sigma,b)
  epsdot=getepsdot(p,sigma,b)
  for i=1:2  
    p[:,i]=rk4!(jefferysRHS,1e-4,2,p[:,i],zeros(3,3),epsdot,1)
    end
  return (p,epsdot)
  end

function getRandc(n)
  c=rand(3,n)-0.5
  for i=1:n
     c[:,i]=c[:,i]/norm(c[:,i])
     end
  return c
  end
getRandc()=getRandc(100)
  c=getRandc()
  epsdot=getepsdot(c,sigma) 
  

function evgroup(c,epsdot,sigma)
  for j=1:10
    epsdot=getepsdot(c,sigma)
    n=size(c)[2]
    for i=1:n
      c[:,i]=rk4!(jefferysRHS,1e-3,n,c[:,i],zeros(3,3),epsdot,1)
      end
    end
  end

function getfps(sigma)
  fp=Array(Float64,3,100)
  for i=1:100
    c=getRandc()
    epsdot=getepsdot(c,sigma)
    evgroup(c,epsdot,sigma)
    fp[:,i]=c[:,1]
    print(i)
    end
  return fp
  end   

function rotC(R)
  # form the K matrix (based on Bowers 'Applied Mechanics of Solids', Chapter 3)
  K1 = [ R[1,1]^2 R[1,2]^2 R[1,3]^2 ; 
         R[2,1]^2 R[2,2]^2 R[2,3]^2 ; 
         R[3,1]^2 R[3,2]^2 R[3,3]^2 ] ; 
  K2 = [ R[1,2]*R[1,3] R[1,3]*R[1,1] R[1,1]*R[1,2] ; 
         R[2,2]*R[2,3] R[2,3]*R[2,1] R[2,1]*R[2,2] ;
         R[3,2]*R[3,3] R[3,3]*R[3,1] R[3,1]*R[3,2] ] ; 
  K3 = [ R[2,1]*R[3,1] R[2,2]*R[3,2] R[2,3]*R[3,3] ; 
         R[3,1]*R[1,1] R[3,2]*R[1,2] R[3,3]*R[1,3] ; 
         R[1,1]*R[2,1] R[1,2]*R[2,2] R[1,3]*R[2,3] ] ; 
  K4 = [ R[2,2]*R[3,3]+R[2,3]*R[3,2] R[2,3]*R[3,1]+R[2,1]*R[3,3] R[2,1]*R[3,2]+R[2,2]*R[3,1] ;
         R[3,2]*R[1,3]+R[3,3]*R[1,2] R[3,3]*R[1,1]+R[3,1]*R[1,3] R[3,1]*R[1,2]+R[3,2]*R[1,1] ;   
         R[1,2].*R[2,3]+R[1,3].*R[2,2] R[1,3].*R[2,1]+ R[1,1].*R[2,3] R[1,1].*R[2,2]+R[1,2].*R[2,1]] ; 
  K = [ K1  2*K2 ; 
        K3   K4   ] ; 
  C = eye(6)
  C[5,5]=10. 
  C[4,4]=10.
  C = K * C * transpose(K) 
  end 

function rotC2(R)
   # form the K matrix (based on Bowers 'Applied Mechanics of Solids', Chapter 3)
  K1 = [ R[1,1]^2 R[1,2]^2 R[1,3]^2 ; 
         R[2,1]^2 R[2,2]^2 R[2,3]^2 ; 
         R[3,1]^2 R[3,2]^2 R[3,3]^2 ] ; 
  K2 = [ R[1,2]*R[1,3] R[1,3]*R[1,1] R[1,1]*R[1,2] ; 
         R[2,2]*R[2,3] R[2,3]*R[2,1] R[2,1]*R[2,2] ;
         R[3,2]*R[3,3] R[3,3]*R[3,1] R[3,1]*R[3,2] ] ; 
  K3 = [ R[2,1]*R[3,1] R[2,2]*R[3,2] R[2,3]*R[3,3] ; 
         R[3,1]*R[1,1] R[3,2]*R[1,2] R[3,3]*R[1,3] ; 
         R[1,1]*R[2,1] R[1,2]*R[2,2] R[1,3]*R[2,3] ] ; 
  K4 = [ R[2,2]*R[3,3]+R[2,3]*R[3,2] R[2,3]*R[3,1]+R[2,1]*R[3,3] R[2,1]*R[3,2]+R[2,2]*R[3,1] ;
         R[3,2]*R[1,3]+R[3,3]*R[1,2] R[3,3]*R[1,1]+R[3,1]*R[1,3] R[3,1]*R[1,2]+R[3,2]*R[1,1] ;   
         R[1,2].*R[2,3]+R[1,3].*R[2,2] R[1,3].*R[2,1]+ R[1,1].*R[2,3] R[1,1].*R[2,2]+R[1,2].*R[2,1]] ; 
  K = [ K1  2*K2 ; 
        K3   K4   ] ; 
  C = zeros(6,6)
  C[5,5]=1. 
  C[4,4]=1.
  C = K * C * transpose(K) 
 

function jefferysRHS(c,vort,epsdot)
  return (vort*c + epsdot*c-(c'*epsdot*c)[1]*c)
  end 

function getRotM(x)
  (norm(x-[0,0,1])<1e-6)?(return eye(3)):pass
  cpMat(v)=[0 v[3] -v[2];-v[3] 0 v[1]; v[2] -v[1] 0]
  z=[0,0,1]
  v=cpMat(x)*z
  s=norm(v)
  c=x[3]
  V=cpMat(v)
  V[1,2]=V[1,2][1]
  V[1,3]=V[1,3][1]
  V[2,3]=V[2,3][1]
  V[2,1]=V[2,1][1]
  V[3,1]=V[3,1][1]
  V[3,2]=V[3,2][1]
   
  return eye(3)+V+(V*V)*(1-c)/s^2  
  end

getRotMinv(x)=inv(getRotM(x))

function getStrW(p,sigma)
  n=size(p)[2]
  Wg=zeros(3,3)
  Dgv=zeros(6)
  D=zeros(3,3)
  Dg=zeros(3,3)
  W=zeros(3,3)
  sigmag=zeros(3,3)
  for i=1:n 
    R=getRotM(p[:,i])
    sigmag=tensor2Voigt(R'*(voigt2Tensor(sigma))*R)
    Dgv[5]=sigmag[5];Dgv[4]=sigmag[4]
    Wg[2,3]=Dgv[4]
    Wg[3,2]=-Dgv[4]
    Wg[1,3]=Dgv[5]
    Wg[3,1]=-Dgv[5]
    Dg=voigt2Tensor(Dgv)
    W+=R*Wg*R'
    D+=R*Dg*R'
    end
  D/=n;W/=n
  return (W,D)
  end

######
p=getRandc()

function getfps(sigma)
  fp=Array(Float64,3,10)
  nt=100
  cr=zeros(3,100)
  for i=1:10
    c=getRandc()
    for j=1:nt
      (W,epsdot)=getStrW(c,sigma)
      evgroup(c,epsdot,sigma,W)
      fp[:,i]=c[:,i]
      cr=c
      end
    print(i)
    end
  return (fp,cr)
  end   

function evgroup(c,epsdot,sigma,W)
    n=size(c)[2]
    for i=1:n
      c[:,i]=rk4!(jefferysRHS,1e-3,n,c[:,i],W,epsdot)
      end
    end

function getepsdot(p,sigma,b)
  C0=rotC(getRotM(p[:,1]))
  C1=rotC(getRotM(p[:,2]))
  C=(1-b)*C0+b*C1
  eps_voigt=C*sigma
  print(eps_voigt)
  return voigt2Tensor(eps_voigt)
  end

function getepsdot(p,sigma)
  n=size(p)[2]
  C=zeros(6,6)
  for i=1:n
    C+=rotC(getRotM(p[:,1]))/n
    end
  eps_voigt=C*sigma
  return voigt2Tensor(eps_voigt)
  end

function getepsW(p,sigma)
  n=size(p)[2]
  C=zeros(6,6)
  for i=1:n
    C+=rotC(getRotM(p[:,1]))/n
    end
    
function voigt2Tensor(v)
  x=zeros(3,3)
  x[1,1]=v[1];x[1,2]=v[6];x[1,3]=v[5]
  x[2,3]=v[4];x[2,2]=v[2];x[3,3]=v[3]
  return symmetrize!(x)
  end
function tensor2Voigt(v)
  x=zeros(6)
  x[1]=v[1,1];x[2]=v[2,2];x[3]=v[3,3]
  x[4]=v[2,3];x[5]=v[1,3];x[6]=v[1,2]
  return x
  end
