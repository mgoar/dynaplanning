module DynaUtils

include("DynaPlanning.jl")

using LightGraphs, LightGraphsFlows, SimpleWeightedGraphs
using GraphPlot, Plots
using Colors
using Compose, Cairo
using SparseArrays
using DataStructures
using Random
using Logging

import .DynaPlanning

"""
price_of_anarchy(r, c)

Compute price of anarchy ρ(G,r,c) for instances with linear cost functions. r and c are vectors of traffic rates and cost functions, respectively.

"""
function price_of_anarchy(r, c)

    C_f = sum([c[i](r[i]) * r[i] for i = 1:length(r)])

    # Optimal flow for the instance (G,r/2,c)
    r_star = 0.5 * r
    C_f_star = sum([c[i](r_star[i]) * r_star[i] for i = 1:length(r_star)])

    # Augmentation to (G,r,c)
    δ = 1
    C_f_star_augm = δ * (C_f - C_f_star)

    return C_f / C_f_star_augm

end

"""
    find_flow(G)

Compute the edge flow for graph `G`.

"""
function find_flow(SInstance)

    A = adjacency_matrix(SInstance.network)

    I, J, vals = findnz(A)

    # init f_e
    f_e = sparse(I, J, zeros(length(vals)))

    for a in SInstance.agents
        thisAgent = a
        thisPath = thisAgent.path
        edge_k = first(thisPath)
        f_e[edge_k.src, edge_k.dst] += 1.0
    end

    # normalize
    f_e = sparse(I, J, findnz(A)[3] / sum(findnz(A)[3]))

    return f_e
end

"""
    find_supply_demand(G)
    
Compute the supply and the demand for `SInstance`.
   
"""
function find_supply_demand(SInstance)

    # init d
    d = zeros(nv(SInstance.network), 1)

    for a in SInstance.agents
        thisAgent = a
        thisPath = thisAgent.path

        @debug thisPath

        edge_src = first(thisPath)
        d[edge_src.src] += 1.0

        edge_sink = last(thisPath)
        d[edge_sink.dst] -= 1.0
    end

    # normalize
    total_demand = sum([x] for x in d if x > zero(x))

    d *= 1 / total_demand

    @debug "Supply/demand:" d

    return d
end

"""
    recursive_DFS!(genericgraph, visited, sink)

Compute -recursively- the list of paths between `visited` (source) and `sink` (destination) pairs for `genericgraph`.

"""
function recursive_DFS!(pathList, genericgraph, visited, sink)

    nodes = LightGraphs.neighbors(genericgraph, visited[length(visited)])

    for i in nodes
        if i == sink
            push!(visited, i)
            push!(pathList, visited)
            pop!(visited)
        elseif !in(i, visited)
            push!(visited, i)
            recursive_DFS!(pathList, genericgraph, visited, sink)
            pop!(visited)
        end
    end

end

"""
    find_all_paths(G, src, dst)

Compute -recursively- all paths between `src` (source) and `dst` (destination) pairs for graph `G`.

"""
function find_all_paths(G, src, dst)

    A = adjacency_matrix(G)

    all_paths = Array{Array{Any},1}()

    DynaPlanning.recursive_DFS!(all_paths, SimpleDiGraph(A), [src], dst)

    return all_paths

end

function plot_graph(G, locs_x, locs_y)

    A = adjacency_matrix(G)

    nodelabel = collect(1:nv(SimpleDiGraph(A)))

    draw(
        PNG("graph.png", 10cm, 10cm),
        gplot(SimpleDiGraph(A), locs_x, locs_y, nodelabel = nodelabel),
    )

end

end