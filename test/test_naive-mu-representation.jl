using PhyloNetworks, Test
include(joinpath(@__DIR__, "..", "mu-representation.jl"))

tre1 = readTopology("((A, B), (C, D));")
map1 = Dict([leaf.name => leaf for leaf in tre1.leaf])
mu1 = node_μ_naive(tre1)

@test mu1[map1["A"]] == 3
@test mu1[map1["B"]] == 3
@test mu1[map1["C"]] == 3
@test mu1[map1["C"]] == 3


net2 = readTopology("((A)#H1,#H1)p;")
map2 = Dict([leaf.name => leaf for leaf in net2.leaf])
mu2 = node_μ_naive(net2)

@test mu2[map2["A"]] == 4


net3 = readTopology("(((A,#H2),#H1),((((B)#H2,#H3))#H1,((#H4)#H3,((C)#H4,D))));")
map3 = Dict([leaf.name => leaf for leaf in net3.leaf])
mu3 = node_μ_naive(net3)

@test mu3[map3["A"]] == 4
@test mu3[map3["B"]] == 11
@test mu3[map3["C"]] == 16
@test mu3[map3["D"]] == 5