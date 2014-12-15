using Compose

function schmidt(p::Array{Float64,2})
    xs=sqrt(2./(1-p'[:,3])).*p'[:,1]
    ys=sqrt(2./(1-p'[:,3])).*p'[:,2]
    return sp_compose((xs/2+1)/2,(ys/2+1)/2)
end

#function sp_compose(xs::Array{Float64,1}, ys::Array{Float64,1})
#    return compose(context(),(context(),circle(xs,ys,0.005*ones(length(xs)))),
#                   (context(0.5w,0.5h),circle(0,0,0.5),fill("bisque"),linewidth(0.1mm),stroke("black")))
#end

function sp_compose(xs::Array{Float64,1}, ys::Array{Float64,1})
    return compose(context(),(context(0,0,1,1),circle(xs,ys,0.005*ones(length(xs)))),
                   (context(0,0,1,1),circle(),fill("bisque"),linewidth(0.1mm),stroke("black")))
end



function squarest(n)
    small=floor(sqrt(n))
    while small>=1
        big=n/small
        if big-int(big)==0
            return (int(big),int(small))
        end
        small-=1
    end
end

function call_by_arr(fun, arr)
    append!(fun.args,arr)
    return eval(fun)
end


function compose_subplots(subplots,m, width)
    m=length(subplots)
    big,small = squarest(m)
    height=small/big*width
    centers_x = rep_els(0:big-1,small)/big
    centers_y = repmat(0:small-1,big)/small
    scale_x = ones(m)/big
    scale_y = scale_x 
    context_arr=Array(Any,m)
    for i=1:m
        context_arr[i] = 
            compose(context(centers_x[i], centers_y[i],
                            scale_x[i], scale_y[i]),
                    subplots[i])
    end
    #return compose(context(0,0,1,1),context_arr)
    return call_by_arr(:(compose(context())),context_arr)
end
