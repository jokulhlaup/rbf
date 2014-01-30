cd("C-axisdatabase")
dr=readdir()
c=Dict()
for d=1:length(dr)
  cd(dr[d])
  rf=readdlm("c-axes.txt",'\n')
  c[int(dr[d])]=readif(rf)
  cd("..")
  end

function readif(rf)
  c=Array(Float64,0)
  for i=1:length(rf)
    
    try  append!(c,float64(split(replace(rf[500],r"[\t\r]","\r"))[2:3]))  #append!(c,float64(split(rf[i],'\t')[2:3]))
      catch
        continue
      end
    end
  return reshape(c,(2,int(length(c)/2)))
  end
