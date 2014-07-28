import Base.indmin
include("SetUtils.jl")

function apAnneal(cost_mat::Array{Float64,2},maxiter::Int)
    n=size(cost_mat,2)
    cv=reshape(deepcopy(cost_mat),n*n)
    (c,w)=greedy!(cv)
    best=w
    temp=2.
    lowtemp=0.0000001
    
    for i=1:maxiter
        temp=1/(3*log(i+1))
        tryProposal!(w,n,cost_mat,temp)
    end
    return w,getCost(cost_mat,best)
end


function getCost(cost_mat,w)
    n=length(w)
    tc=0.
    for i=1:n
        tc+=cost_mat[w[i],i]
    end
    return tc
end
function tryProposal!(w::Array{Int,1},
                      n::Int,
                      cost_mat::Array{Float64,2},
                      temp::Float64,
                      )
    j=rand(1:n);k=rand(1:n)
    w2=deepcopy(w)
    w2[k]=w[j];w2[j]=w[k]
#    dcost=getCost(cost_mat,w2)-getCost(cost_mat,w)
    dcost=cost_mat[w[k],j]+cost_mat[w[j],k]
          -cost_mat[w[j],j]-cost_mat[w[k],k]
    if accept!(temp,dcost)==true
        rev2!(w,n,j,k)
    end
end

function sort!(v::AbstractVector,inds, lo::Int, hi::Int)
    SMALL_THRESHOLD=1e-20
    @inbounds while lo < hi
        hi-lo <= SMALL_THRESHOLD && return sort!(v,inds, lo, hi)
        mi = (lo+hi)>>>1
        if <( v[mi], v[lo])
            v[mi], v[lo] = v[lo], v[mi]
            inds[mi], inds[lo] = inds[lo], inds[mi]
        end
        if <( v[hi], v[mi])
            v[hi], v[mi] = v[mi], v[hi]
            inds[hi], inds[mi] = inds[mi], inds[hi]
        end
        if <( v[mi], v[lo])
            v[mi], v[lo] = v[lo], v[mi]
            inds[mi], inds[lo] = inds[lo], inds[mi]
        end
        v[mi], v[lo] = v[lo], v[mi]
        inds[mi], inds[lo] = inds[lo], inds[mi]
        i, j = lo, hi
        pivot = v[lo]
        while true
            i += 1; j -= 1;
            while <( v[i], pivot); i += 1; end;
            while <( pivot, v[j]); j -= 1; end;
            i >= j && break
            v[i], v[j] = v[j], v[i]
            inds[i], inds[j] = inds[j], inds[i]
        end
        v[j], v[lo] = v[lo], v[j]
        inds[j], inds[lo] = inds[lo], inds[j]
        lo < (j-1) && sort!(v,inds, lo, j-1)
        lo = j+1
    end
    return v
end




function accept!(temp,dcost)
    if dcost<0
        return true
    elseif (exp(-dcost/temp) > rand())
        return true
    else
        return false
    end
end 

flat2mat(i,n)=int((i-1)%n+1),int(ceil(i/n))


function greedy!(cost::Vector)
    n=length(cost)
    m=int(sqrt(n))
    res=Array(Float64,m)
    inds=int(linspace(1,n,n))
    sort!(cost,inds,1,n)
    w=Array(Int,m)
    avail_cols=IntSet(1:m)
    avail_inds=IntSet(1:m)
    navail=n
    counter=n
    tc=0.
    i=1
    while counter > 0
        (l,k)=flat2mat(inds[i],m) 
        if ((k in avail_cols) & (l in avail_inds))
            res[k]=l
            tc+=cost[i]
            delete!(avail_cols,k)
            delete!(avail_inds,l)
        end
        counter -= 1
        i += 1
        end
    return tc,int(res)
end


    



function transition!(w::Array{Int,1},n::Int)
    i=rand(1:n);j=rand(1:n)
    rev2!(w,n,i,j)
    return (i,j)
end

function rev2!(w::Array{Int,1},
              n::Int,
              i::Int,
              j::Int)
    old1=w[i]
    old2=w[j]
    w[j]=old1
    w[i]=old2

end

function getInd(w::Array{Int,1},ind,n)
    n=length(w)
    for i=1:n
        if w[i]==ind
            return i
        end
    end
    error("No such index $ind in w")
end


function rand_xs(seed::Int64)
    seed >>= 12
    seed <<= 25
    seed >>= 27
    return seed * 2685821657736338717
end
