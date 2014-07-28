using PyCall
@pyimport hungarian
include("SA.jl")

function assignment(C)
    inds=hungarian.lap(C)[1] .+ 1
    tc=getCost(C,inds)
    return (tc,inds)
end

function assignment!(C,S,n)
    @assert typeof(S)==Array{Int32,1}
    @assert size(S)==(n,)
    @assert size(C)==(n*n,)
    @assert typeof(C)==Array{Float64,1}

    mincost=ccall((:solve_LSAP,"assignment"), Cdouble, (Ptr{Cdouble},Cint,Ptr{Cint}), C,n,S)
    return mincost
end

#t = ccall( (:solve_la, "assignment"), Int32,(Ptr{Ptr{Cdouble}},Int32,Ptr{Cdouble},Ptr{Cint}),&C,n,&cost,S)
 
