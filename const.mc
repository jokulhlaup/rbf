cpMat(v):=-matrix([0,v[3],-v[2]],[-v[3],0,v[1]],[v[2],-v[1],0]);
norm(x):=sqrt(x.x);

/*Rodrigues rotation formula*/
R(x):=(
    z:[0,0,1],
    v:cpMat(x).z,
    s:norm(v),
    c:x[3],
    V:cpMat(v),
    V[1,2]:V[1,2][1],
    V[1,3]:V[1,3][1],
    V[2,3]:V[2,3][1],
    V[2,1]:V[2,1][1],
    V[3,1]:V[3,1][1],
    V[3,2]:V[3,2][1],
   
    ident(3)+V+(V.V)*(1-c)/s^2)

/*rotate the fluidity tensor by R in Voigt notation*/
rotC(R):=(   K1 : transpose(mat_unblocker(matrix(
                 [ R[1,1]^2,R[1,2]^2, R[1,3]^2 ], 
                 [ R[2,1]^2, R[2,2]^2, R[2,3]^2 ], 
                 [ R[3,1]^2, R[3,2]^2, R[3,3]^2 ] ))),

             K2 : transpose(mat_unblocker(matrix(
                 [ R[1,2]*R[1,3], R[1,3]*R[1,1], R[1,1]*R[1,2] ],
                 [ R[2,2]*R[2,3], R[2,3]*R[2,1], R[2,1]*R[2,2] ],
                 [ R[3,2]*R[3,3], R[3,3]*R[3,1], R[3,1]*R[3,2] ] ))),

             K3: transpose(mat_unblocker(matrix(
                 [ R[2,1]*R[3,1], R[2,2]*R[3,2], R[2,3]*R[3,3] ],
                 [ R[3,1]*R[1,1], R[3,2]*R[1,2], R[3,3]*R[1,3] ], 
                 [ R[1,1]*R[2,1], R[1,2]*R[2,2], R[1,3]*R[2,3] ] ))),

             K4: transpose(mat_unblocker(matrix(
                  [ R[2,2]*R[3,3]+R[2,3]*R[3,2], R[2,3]*R[3,1]+R[2,1]*R[3,3], 
                    R[2,1]*R[3,2]+R[2,2]*R[3,1] ], 
                  [  R[3,2]*R[1,3]+R[3,3]*R[1,2], R[3,3]*R[1,1]+R[3,1]*R[1,3], 
                    R[3,1]*R[1,2]+R[3,2]*R[1,1] ],    
                  [  R[1,2]*R[2,3]+R[1,3]*R[2,2], R[1,3]*R[2,1]+
                    R[1,1]*R[2,3], R[1,1]*R[2,2]+R[1,2]*R[2,1] ] ))),

             K:transpose(mat_unblocker(matrix(
                 [K1, K3],[2*K2,K4]))),
             C:ident(6),
             C[5,5]:hf,
             C[4,4]:hf,
             K.C.transpose(K))

/* Maxima is the worst. */
dDot(a,b):= (ab: a*b,
             sum(transpose(sum(row(ab,i),i,1,3))[j],j,1,3)[1]
             )
outer(a,b):=matrix(a[1]*b,a[2]*b,a[3]*b)

jefferysRHS(D,W,p):= W.p+D*p-dDot(D,outer(p,p))*p

strainRate(S_voigt,C) :=( str0:C.S_voigt,
                          s:str0*[1,1,1,1/2,1/2,1/2],
                          strR:matrix([s[1],s[6],s[5]],[s[6],s[2],s[4]],[s[5],s[4],s[3]])
                          )
C_p:matrix([1,0,0,0,0,0],[0,1,0,-2*a*f+a,0,0],
           [0,0,1,2*a*f-a,0,0],[0,-2*a*f+a,2*a*f-a,f,0,0],
           [0,0,0,0,f,a*(1-f)],[0,0,0,0,a*(1-f),1]);

C:ident(6);
C[5,5]:hf;
C[4,4]:hf;

#C perturbed by another population of crystals
C_p2:(b*rotC(R([1,q,0]))+C*(1-b));
S:[1,-1/2,-1/2,0,0,0];

[vals,vects]=eigenvectors(Voigt2Tensor(C_p2.S))

#perturbed D in matrix form 
Voigt2Tensor(Sp):=matrix([Sp[1,1],Sp[6,1],Sp[5,1]],[Sp[6,1],Sp[2,1],Sp[4,1]],[Sp[5,1],Sp[6,1],Sp[3,1]]);           
Dp2=Voigt2Tensor(Sp);

dcDp2:diff(rotC(R([1,a,b])),a);
Jp2:ev(dcDp2,a=0,b=0)
dcDp1:diff(rotC(R([1,a,b])),b);
Jp1:ev(dcDp1,a=0,b=0)


unbox(m):=matrix([m[1,1][1],m[1,2][1],m[1,3][1]],[m[2,1][1],m[2,2][1],m[2,3][1]],[m[3,1][1],m[3,2][1],m[3,3][1]])

objective(p):=jefferysRHS(strainRate(S_voigt,rotC(R(p))),0,p))
