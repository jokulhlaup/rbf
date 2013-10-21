using jefferys
ngr=100
ns=10
coors=rand(ns,3)
p=rand(10,100,3)
p=p./(sqrt(p[:,:,1].^2+p[:,:,2].^2+p[:,:,3].^2))
h=1
C=Array(Float64,ns,6,6)
fab=Fabric{Float64,Int64}(coors,p,ngr,ns,h,C)



