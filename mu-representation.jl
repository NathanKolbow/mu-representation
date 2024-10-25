using PhyloNetworks, DataStructures

const Node = PhyloNetworks.Node;
const Edge = PhyloNetworks.Edge;
const Labels = AbstractVector{<:AbstractString};


"""
`N` is treated as semi-directed - i.e. all edges are treated as undirected except
hybrid edges which are treated as directed. Leaf edges also treated as undirected.
"""
function edge_μ_semi_directed_naive(N::HybridNetwork; L::Labels=tipLabels(N))

    # Each edge has 1 entry if it's a hybrid, two if it's not.
    # In the tree-edge case, we need to keep track of which direction both entries correspond to.
    mapping = Dict{
        Tuple{Edge, Tuple{Node, Node}}, # relevant edge AND the direction being tracked
                                        # hybrid edges have 1 unique entry, tree edges have 2
        Tuple{Vector{Int}, Char}        # the mapping associated with that edge + direction
    }();
    label_map = Dict{String, Int}([taxon => i for (i, taxon) in enumerate(L)])
    logged = Dict{Tuple{Edge, Tuple{Node, Node}}, Bool}()

    # Populate the mapping w/ initial values
    for edge in N.edge
        if edge.hybrid
            direction = (getparent(edge), getchild(edge))
            mapping[(edge, direction)] = (zeros(length(L)), 'h')
            logged[(edge, direction)] = false
        else
            direction = (getparent(edge), getchild(edge))
            mapping[(edge, direction)] = (zeros(length(L)), 't')
            logged[(edge, direction)] = false

            direction = (getchild(edge), getparent(edge))
            mapping[(edge, direction)] = (zeros(length(L)), 't')
            logged[(edge, direction)] = false
        end
    end

    # Populate the queue (naively)
    Q = Queue{Tuple{Edge, Tuple{Node, Node}}}()
    for edge in N.edge
        enqueue!(Q, (edge, (getparent(edge), getchild(edge))))
        if !edge.hybrid
            enqueue!(Q, (edge, (getchild(edge), getparent(edge))))
        end
    end


    j = 0;
    while !isempty(Q)
        j += 1;
        if j > 10000 error("Looped $(j) times - exiting.") end

        edge, (from_node, to_node) = dequeue!(Q)
        direction = (from_node, to_node)

        if to_node.leaf
            # Leaf edge
            # @info "(L) $(from_node.name) --> $(to_node.name)"
            mapping[(edge, direction)][1][label_map[to_node.name]] += 1
            logged[(edge, direction)] = true
            continue;
        else
            path_edges = [E for E in to_node.edge if E != edge && (!E.hybrid || getparent(E) == to_node)]
            if any(E -> !logged[E, (to_node, other_node(E, to_node))], path_edges)
                # We don't have enough other edge data to log this edge in this direction yet,
                # so re-queue this edge + direction and continue for now
                enqueue!(Q, (edge, direction))
                continue
            end

            # @info "(T) $(from_node.name) --> $(to_node.name)"
            for E in path_edges
                E_dir = (to_node, other_node(E, to_node))
                mapping[(edge, direction)][1] .+= mapping[(E, E_dir)][1]
            end
            logged[(edge, direction)] = true
        end
    end

    # Still need the \mu_R(T) bit
    return collect(values(mapping))

end
tre1 = readTopology("((A,B)ab,(C,D)cd)r;");
tre1_map = edge_μ_semi_directed_naive(tre1);
edge_μ_dist(tre1, tre1, true)

N1 = readTopology("((a,#H1)aa,b,(((c)cc)#H1,d)cd)r;")
N2 = readTopology("((a,#H1),b,(((c))#H1,d));")
edge_μ_semi_directed_naive(N1)
edge_μ_dist(N1, N2, true)




function other_node(edge::Edge, node::Node)
    return edge.node[1] == node ? edge.node[2] : edge.node[1]
end


# Naive algorithm for computing the node-based μ-representation of network `N`.
#
# Note: treats `N` as ROOTED, i.e. the only "root component" is trivial: the root itself
#
# Note: all edges in `N` are treated as directed! Thus, for network (A,B) we
# would say there are **0** paths from A to B or vice versa.
#
# The node-based μ-representation of `N` is as follows:
#     \$\$\mu_V(N) = \int \mu(v, N); v\in V(N)\$\$
# where \$V(N)\$ is the set of vertices in `N` and \$\mu(v, N)\$ is the tuple \$(\mu_1(v), \dots, \mu_n(v))\$ where
# \$\mu_i(v)\$ is the number of paths in `N` from \$v\$ to the i'th leaf of `N`.
#
# Returns a `Dict{PhyloNetworks.Node, Int}` object mapping each leaf \$l_i\$ in `N` to \$\sum_{v\in V} \mu(v,N)[i]\$,
#    i.e. the number of paths from a given node to leaf \$i\$ summed across all nodes in `N`.
function node_μ_naive(N::HybridNetwork)

    L = tipLabels(N)
    mapping = Dict{PhyloNetworks.Node, Int}([leaf => 0 for leaf in N.leaf])

    for node in N.node
        for desc_leaf in get_descendant_leaves(node, multiplicity=true)
            mapping[desc_leaf] += 1
        end
    end

    return mapping

end


function edge_μ_naive(N::HybridNetwork; L::AbstractVector{<:AbstractString}=tipLabels(N), unrooted::Bool=false)
    # {(μ_r(T), :r)} : T \in R is exactly \mu(r, N), r the root node
    !unrooted || error("NOT IMPLEMENTED FOR UNROOTED TREES/NETWORKS")

    leaf_idx = Dict{AbstractString,Int}([leaf => i for (i, leaf) in enumerate(L)])
    mapping = Dict{Union{Edge, Node}, Tuple{Vector{Int},Char}}([e => (zeros(N.numTaxa), ifelse(e.hybrid, 'h', 't')) for e in N.edge])

    for edge in N.edge
        # \mu(edge) = \mu(getchild(edge), N), which we get w/ `get_descendant_leaves` with multiplicity
        for desc_leaf in get_descendant_leaves(getchild(edge), multiplicity=true)
            mapping[edge][1][leaf_idx[desc_leaf.name]] += 1
        end
    end

    # then, since we assume rooted, for the rooted components, just add the root's \mu(root, N)
    mapping[N.node[N.root]] = (zeros(N.numTaxa), 'r')
    for desc_leaf in get_descendant_leaves(N.node[N.root], multiplicity=true)
        mapping[N.node[N.root]][1][leaf_idx[desc_leaf.name]] += 1
    end

    return collect(values(mapping))
end


function edge_μ_dist(N1::HybridNetwork, N2::HybridNetwork, semi_directed::Bool=true)
    (N1.numTaxa == N2.numTaxa && all(t in tipLabels(N2) for t in tipLabels(N1))) || error("N1 and N2 must be defined on the same leaf set.")

    L = tipLabels(N1)
    N1_eμ = semi_directed ? edge_μ_semi_directed_naive(N1; L=L) : edge_μ_naive(N1; L=L)
    N2_eμ = semi_directed ? edge_μ_semi_directed_naive(N2; L=L) : edge_μ_naive(N2; L=L)

    return length(symdiff(N1_eμ, N2_eμ))
end


"""
`multiplicity`: whether leaves should be returned once for each path that
                leads to them from `node`, or only ever once
"""
function get_descendant_leaves(node::PhyloNetworks.Node; multiplicity::Bool=false)
    q = Vector{PhyloNetworks.Node}([node])
    leaves = Vector{PhyloNetworks.Node}([])
    while length(q) > 0
        curr = pop!(q)
        if curr.leaf
            push!(leaves, curr)
        else
            for child in getchildren(curr)
                push!(q, child)
            end
        end
    end
    
    if multiplicity return leaves end
    return unique(leaves)
end