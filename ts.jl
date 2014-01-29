cd("C-axisdatabase")
dr=readdir()
c=Dict()
for d in dr
  cd(d)
  rf=readdlm("c-axes.txt",'\n')
  c[int(d)]=readif(rf)
  cd("..")
  end

function readif(rf)
  c=Array(Float64,0)
  for i=1:length(rf)
    try append!(c,float64(split(rf[i],'\t')[2:3]))
      catch
        continue
      end
    end
  return reshape(c,(2,int(length(c)/2)))
  end
