cost=rand(100)
x=deepcopy(cost)
tc,w=greedy!(x)
best=w
cost=reshape(cost,10,10)
tryProposal!(w,10,cost,1.,best)

w=apAnneal(cost,10000)

function neighbor!(w,w2)
    w2=rev2!(w,10,rand(1:10),rand(1:10))
end
