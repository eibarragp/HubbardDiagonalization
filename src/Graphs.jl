"Defines a basic undirected graph structure and some common graph types."
module Graphs

export Graph, from_cluster, num_sites, edges, has_edge, linear_chain

"Represents an undirected graph."
struct Graph
    num_sites::Int
    # These edges are stored as (min, max) tuples for easy look-up
    edges::Set{Tuple{Int,Int}}
end

function from_cluster(cluster::Tuple)
	num_sites = length(cluster[1])
	edges = Set{Tuple{Int,Int}}()
	for edge in cluster[2]
		push!(edges, (min(edge[1], edge[2]), max(edge[1], edge[2])))
	end

	return Graphs.Graph(num_sites, edges)
end

"""
	linear_chain(n_sites::Int)
Generates a graph with `n_sites` sites organized into a linear chain with periodic boundary conditions.
"""
function linear_chain(n_sites::Int)
    edges = Set{Tuple{Int,Int}}()
    for i in 1:(n_sites-1)
        push!(edges, (i, i+1))
    end
    if n_sites > 2
        push!(edges, (1, n_sites))  # Periodic boundary conditions
    end
    return Graph(n_sites, edges)
end

"""
	num_sites(g::Graph)
Returns the number of sites in the graph `g`.
"""
function num_sites(g::Graph)
    return g.num_sites
end

"""
	edges(g::Graph)
Returns the set of edges in the graph `g`.
Each edge is represented as a tuple `(i, j)` with `i < j`.
"""
function edges(g::Graph)
    return g.edges
end

"""
	has_edge(g::Graph, i::Int, j::Int)
Returns `true` if there is an edge between sites `i` and `j` in the graph `g`, `false` otherwise.
"""
function has_edge(g::Graph, i::Int, j::Int)
    return (min(i, j), max(i, j)) in g.edges
end

end
