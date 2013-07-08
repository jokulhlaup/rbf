#Generate a sparse distance matrix

def getDistMat(xs,p)
   (d,indices)=xs.query(xs)

   for i in xrange(d.shape[1])
      #
