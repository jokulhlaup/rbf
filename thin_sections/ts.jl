cd("C-axisdatabase")
dr=readdir()
c=Dict()
p=Dict()
svs=Array(Float64,length(dr))
for d=1:length(dr)
  cd(dr[d])
  rf=readdlm("c-axes.txt",'\n')
  i=int(dr[d])
  c[i]=readif(rf)
  #convert to rads
  c[i]=c[i]*pi/180
  p[i]=Array(Float64,3,length(c[i][1,:]))
  p[i][1,:]=cos(c[i][1,:])
  p[i][2,:]=sin(c[i][1,:]).*cos(c[i][2,:])
  p[i][3,:]=sin(c[i][1,:]).*sin(c[i][2,:])
  x=svd(p[i])[2]
  svs[d]=x[1]/norm(x)
  cd("..")
  end

function readif(rf)
  c=Array(Float64,0)
  for i=1:length(rf)
    
    try  append!(c,float64(split(replace(rf[i],r"[\t\r]","\r"))[2:3]))  #append!(c,float64(split(rf[i],'\t')[2:3]))
      catch
        continue
      end
    end
  return reshape(c,(2,int(length(c)/2)))
  end

#get singular values and directions
