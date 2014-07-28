function indmin(s::IntSet,
                row)
   rowm=1e20
   for a in s
       if x[a]<rowm
           rowm=x[a]
       end
   end
   return rowm
end

function nth(s::IntSet,
             n::Int)
    i=1
    for a in s
       if i==n
           return a  
           else i+=1
           end
    end
    error("$n out of bounds")
end
       
