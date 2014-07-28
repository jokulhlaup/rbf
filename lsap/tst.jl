using PyCall
@pyimport hungarian
include("SA.jl")

function assignment(C)
    inds=hungarian.lap(C)[1] .+ 1
    return inds
end

function hungarian_stencil(costs,i,j,n)

