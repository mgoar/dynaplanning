module TestDynaPlanning_Instances

include("../src/DynaPlanning.jl")
include("../src/DynaUtils.jl")

using LightGraphs, LightGraphsFlows, GraphPlot, Compose
using SparseArrays
using DataStructures
using Random
using Test
using Logging

import .DynaPlanning
import .DynaUtils

# Logging
logger = LogLevel(Logging.Debug)

mutable struct SimulationParameters
    episodes::Int64
end

parameters = SimulationParameters(100)

################################################################################
# Pigou's example
################################################################################

sources = [1, 2, 1, 4]
sinks = [2, 3, 4, 3]

Pigou = SimpleDiGraph(length(sources))

for i = 1:length(sources)
    add_edge!(Pigou, sources[i], sinks[i])
end

a_e = [0.0, 0.0, 1.0, 0.0]
b_e = [1.0, 0.0, 0.0, 0.0]

# Vector of linear cost functions
c_e = Array{Function}(undef, 4)
for i = 1:length(sources)
    c_e[i] = x -> a_e[i] * x + b_e[i]
end

################################################################################
# # Create SimpleWeightedDiGraph from a Barabási–Albert model random graph

G = barabasi_albert(12, 2, 2, is_directed=true, seed=23)
A = adjacency_matrix(G)

# Cost function
rng = MersenneTwister(1234);
a_e_G = rand!(rng, zeros(12))

rng = MersenneTwister(5678);
b_e_G = rand!(rng, zeros(12))

locs_x, locs_y = spring_layout(G)

@info("Plotting Barabási–Albert network (G)")
DynaUtils.plot_graph(G, locs_x, locs_y)

# DynaAgents
N_DynaAgents = 4

@info "Setting Dyna_agentPlus and StackelbergInstance"  N_DynaAgents 

srcDst = Array{Tuple,1}()
push!(srcDst, (11, 1))
push!(srcDst, (12, 1))
push!(srcDst, (12, 1))
push!(srcDst, (10, 1))

agents = Array{DynaPlanning.DynaAgentPlus,1}()
for n in range(1, stop=N_DynaAgents)
    testAgent = DynaPlanning.DynaAgentPlus()
    testSrcDst = pop!(srcDst)

    # A* shortest path
    testShortestPath = a_star(SimpleDiGraph(A), testSrcDst[1], testSrcDst[2])
    
    testAgent.path = Queue{Edge}()
    for i in 1:length(testShortestPath)
        enqueue!(testAgent.path, testShortestPath[i])
    end

    push!(agents, testAgent)
end

# Stackelberg instance
testStackelbergInstance = DynaPlanning.StackelbergInstance()
testStackelbergInstance.network = G
testStackelbergInstance.agents = agents

end # module