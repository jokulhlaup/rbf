module jefferys
export solveJefferys,ode4u
#Modification of ODE4 from package ODE
function ode4u{T}(F::Function, h::Float64,n::Int64, x0::AbstractVector{T})
    x = Array(T, (length(tspan), length(x0)))
    x[1,:] = x0'

    midxdot = Array(T, (4, length(x0)))
    for i = 1:n
        # Compute midstep derivatives
        midxdot[1,:] = F(x[i,:]')
        midxdot[2,:] = F(x[i,:]' + midxdot[1,:]'*h/2)
        midxdot[3,:] = F(x[i,:]' + midxdot[2,:]'*h/2)
        midxdot[4,:] = F(x[i,:]' + midxdot[3,:]'*h)

        # Integrate
        x[i+1,:] = x[i,:] + 1./6.*h[i].*[1 2 2 1]*midxdot
    end
    return x
end


#Using ODE
function solveJefferys{T}(vort,epsdoy,theta,phi,dt,n)
  n=theta->[sin(theta[1])*cos(theta[2]),sin(theta[1])*sin(theta[2]),cos(theta[1])]

  f=theta->((vort*n(theta))[2:3])-((epsdot*n(theta))[2:3]-(n(theta)'*epsdot*n(theta))[0]*n(theta)[2:3])
  n=jefferys.ode4u(f,dt/4,4,theta)
  theta[1]=acos(lhs[3])
  theta[2]=asin(lhs/(asin(theta)))
  return theta
end

end

function foo(theta:)
((vort*n(theta))[2:3])-((epsdot*n(theta))[2:3]-(n(theta)'*epsdot*n(theta))*n(theta)[2:3]) 
end



