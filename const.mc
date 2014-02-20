cpMat(v):=-matrix([0,v[3],-v[2]],[-v[3],0,v[1]],[v[2],-v[1],0])
norm(x):=sqrt(x.x)
R(x):=(V:cpMat(x),
    z:[0,0,1],
    v:cpMat(x).z,
    s:norm(v),
    c:x[3],
    V:cpMat(v),
    ident(3)+V+(V.V)*(1-c)/s^2
    )

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
             C:ident(6)*10^3,
             C[5,5]:1,
             C[4,4]:1,
             K.C.transpose(K))

(defun outer ($x $y)




