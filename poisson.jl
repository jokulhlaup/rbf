function rhs(x,x0,ep,n)
  for i=1:n
    L[:,i]=d2imq(x[:,i],x0[:,i],ep)
    end
  end

function genrStencil
