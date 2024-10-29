using PhyloNetworks, PhyloCoalSimulations, Test
include(joinpath(@__DIR__, "rf.jl"))
include(joinpath(@__DIR__, "..", "mu-representation.jl"))

const NSIM = 50;
star_tre = readTopology("(((((((((((A,B),C),D),E),F),G),H),I),J),K),L);");
for edge in star_tre.edge edge.length = 0.0 end
sim_tres = simulatecoalescent(star_tre, NSIM+1, 1)

for j = 1:NSIM
    @test edge_μ_dist(sim_tres[j], sim_tres[j+1], true) == 2*hardwiredClusterDistance(sim_tres[j], sim_tres[j+1], false)
end

