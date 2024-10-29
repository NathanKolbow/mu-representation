using PhyloNetworks, DataStructures

const Node = PhyloNetworks.Node;
const Edge = PhyloNetworks.Edge;


"""
`N` is treated as semi-directed - i.e. all edges are treated as undirected except
hybrid edges which are treated as directed. Leaf edges also treated as undirected.
"""
function edge_μ_semi_directed(N::HybridNetwork; L::AbstractVector{<:AbstractString}=tipLabels(N))

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
        if j > 100000 error("Looped $(j) times - exiting.") end

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

    # @warn "Mappings for root components are not yet computed"
    return collect(values(mapping))


    # Get the root components of the network, then calculate \mu_R(T) for each root component T
    root_components = gather_root_components(N)
    root_mappings = []
    for T in root_components
        # T <: Vector{Union{Node, Edge}}
        # if trivial, T is the root node
        # if non-trivial, T will have at least 1 edge E = uv, in which case
        # \mu_R(T) is equivalent regardless of which E we choose, so we just pick the first E we see
        if length(T) == 1
            # trivial - T = [<root node>]
            # QUESTION: from my understanding of Lemma 18 and Definition 10, when T is trivial (a single node), \mu_r(T)
            #           is simply \mu(T[1], N) as defined at the beginning of Section IV, i.e. {#1, #2, ..., #n} where
            #           #i is the number of paths connecting the root to leaf i - is this correct?
            typeof(T[1]) <: Node || error("SANITY CHECK FAILED: T[1] has type $(typeof(T[1])) when T has size $(length(T))")
            root_node = T[1]
            root_node_mapping = zeros(length(L))
            for E in root_node.edge
                direction = (root_node, other_node(E, root_node))
                root_node_mapping += mappings[(E, direction)][1]
            end
            push!(root_mappings, (root_node_mapping, 'r'))
        else
            # non-trivial - take the first edge in the set
            any(typeof(comp) <: Edge for comp in T) || error("SANITY CHECK FAILED: no edges in root component of size $(length(T)) ??")
            E = T[findfirst(component -> typeof(component) <: Edge, T)]

            E_mapping = zeros(length(L))
            for direction in ((getparent(E), getchild(E)), (getchild(E), getparent(E)))
                E_mapping += mappings[E, direction][1]
            end
            push!(root_mappings, (E_mapping, 'r'))
        end
    end
    @show root_mappings


    return [collect(values(mapping)); root_mappings]
end


function other_node(edge::Edge, node::Node)
    return edge.node[1] == node ? edge.node[2] : edge.node[1]
end


function gather_root_components(N::HybridNetwork)
    # A root component is an undirected, maximal component of the network
    #
    # At the moment, I am not entirely clear on what is and is not a root component, eg:
    # QUESTION: the leaf node h and its parent nodes are undirected components,
    #           so why are they not root components?
end


function edge_μ_dist(N1::HybridNetwork, N2::HybridNetwork; semi_directed::Bool=true)
    (N1.numTaxa == N2.numTaxa && all(t in tipLabels(N2) for t in tipLabels(N1))) || error("N1 and N2 must be defined on the same leaf set.")

    L = tipLabels(N1)
    N1_eμ = semi_directed ? edge_μ_semi_directed(N1; L=L) : error("only semi-directed allowed")
    N2_eμ = semi_directed ? edge_μ_semi_directed(N2; L=L) : error("only semi-directed allowed")

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