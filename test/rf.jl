using PhyloNetworks


function RobinsonFoulds(tree1::HybridNetwork, tree2::HybridNetwork)
    (tree1.numHybrids == 0 && tree2.numHybrids == 0) || error("tree1 and tree2 must be tree-like (0 hybrids)")
    
    clusters1 = getHardwiredClusters(tree1)
    clusters2 = getHardwiredClusters(tree2)
    inboth = intersect(clusters1, clusters2)

    return length(clusters1) + length(clusters2) - 2 * length(inboth)
end


function getHardwiredClusters(tree::HybridNetwork)
    working_edgeset = [leaf.edge[1] for leaf in tree.leaf]
    saved_map = Dict()      # final clusters are the values of the dict
    rootnumber = tree.node[tree.root].number 

    while length(working_edgeset) != 0
        edge = working_edgeset[1]
        deleteat!(working_edgeset, 1)
        if edge.number in keys(saved_map) continue end
        # @printf "Working edge %d\n" edge.number

        if getchild(edge).leaf > 0
            # Leaf edge
            saved_map[edge.number] = [getchild(edge).name]
            if getparent(edge) != tree.node[tree.root] push!(working_edgeset, getparentedge(getparent(edge))) end
        else
            # Non-leaf edge
            children = getchildren(getchild(edge))
            e1 = getparentedge(children[1])
            e2 = getparentedge(children[2])

            try
                saved_map[edge.number] = vcat(saved_map[e1.number], saved_map[e2.number])
            catch e
                push!(working_edgeset, edge)
            end
            # saved_map[edge.number] = vcat(saved_map[edge.node[2].edge[1].number], saved_map[edge.node[2].edge[2].number])
        
            # Keep going up
            if getparent(edge) != tree.node[tree.root] push!(working_edgeset, getparentedge(getparent(edge))) end
        end

    end

    for key in keys(saved_map) saved_map[key] = sort(saved_map[key]) end
    return values(saved_map)
end