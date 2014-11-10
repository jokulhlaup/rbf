function rhs(x,x0,ep,n)
  for i=1:n
    L[:,i]=(d2imq(x[:,i],x0[:,i],ep,1,1)
      +d2imq(x[:,i],x0[:,i],ep,2,2))
    end
  end

function genrStencil
  #get some points
