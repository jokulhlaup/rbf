#Modification of ODE4 from package ODE
function ode4u{T}(F::Function, h::Float64,n::Int64, x0::AbstractVector{T})
    x = Array(T, (length(tspan), length(x0)))
    x[1,:] = x0'

    midxdot = Array(T, (4, length(x0)))
    for i = 1:n
        # Compute midstep derivatives
        midxdot[1,:] = F(tspan[i], x[i,:]')
        midxdot[2,:] = F(tspan[i]+h/2, x[i,:]' + midxdot[1,:]'*h/2)
        midxdot[3,:] = F(tspan[i]+h/2, x[i,:]' + midxdot[2,:]'*h/2)
        midxdot[4,:] = F(tspan[i]+h, x[i,:]' + midxdot[3,:]'*h)

        # Integrate
        x[i+1,:] = x[i,:] + 1./6.*h[i].*[1 2 2 1]*midxdot
    end
    return (tspan, x)
end


#Using ODE
function solveJefferys(vort::Array{T,2},epsdot:Array{T,2},theta,phi,dt,n)
  n=(theta,phi)->[sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)]

  f=(theta,phi)->vort[2:3,:]*n(theta,phi)[2:3]-(epsdot[2:3,:]*n(theta,phi)[2:3]-(n(theta,phi)*epsdot*n)*n(theta,phi) [2:3])
  





