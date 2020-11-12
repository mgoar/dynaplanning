module testInstances

include("DynaPlanning.jl")

using LightGraphs, LightGraphsFlows, SimpleWeightedGraphs, GraphPlot, Compose # to reduce scope
using SparseArrays                      # idem
using DataStructures                    # idem
using Random
using Test
using Logging

import .DynaPlanning

global rng_1 = MersenneTwister(1989) # seed
global rng_2 = MersenneTwister(2302) # seed

# create SimpleWeightedDiGraph
sources = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 10, 4, 7, 10, 5, 2, 2, 1, 9]
destinations = [2, 8, 5, 5, 1, 10, 2, 4, 5, 7, 2, 8, 10, 6, 9, 9, 1, 10, 2]
wghts = [
    1,
    0.2,
    0.1,
    0.45,
    0.05,
    0.1,
    0.45,
    0.8,
    0.1,
    0.45,
    0.95,
    0.6,
    0.1,
    0.05,
    0.6,
    0.3,
    0.1,
    0.65,
    0.035,
]

G = SimpleWeightedDiGraph(sources, destinations, wghts)
A = adjacency_matrix(G)

locs_x, locs_y = spring_layout(G)

DynaPlanning._util_plot_graph(G, locs_x, locs_y)

testAgent1 = DynaPlanning.Dyna_agentPlus()
testAgent2 = DynaPlanning.Dyna_agentPlus()

testPath1 = a_star(SimpleDiGraph(A), 3, 10, weights(G))
testAgent1.path = Queue{Edge}()
for i in 1:length(testPath1)
     enqueue!(testAgent1.path, testPath1[i])
 end

testPath2 = a_star(SimpleDiGraph(A), 6, 5, weights(G))
testAgent2.path = Queue{Edge}()
 for i in 1:length(testPath2)
      enqueue!(testAgent2.path, testPath2[i])
  end

testStackelbergInstance = DynaPlanning.StackelbergInstance()

testStackelbergInstance.network = G
testStackelbergInstance.β = 0.50
testStackelbergInstance.agents = Array{DynaPlanning.Dyna_agentPlus,1}()
push!(testStackelbergInstance.agents, testAgent1)
push!(testStackelbergInstance.agents, testAgent2)

# linear cost functions coefficients - constant terms
M = LightGraphs.ne(G)
W = LightGraphs.weights(Main.testInstances.testStackelbergInstance.network)
I = findnz(W.parent)[1]
J = findnz(W.parent)[2]

b_i = sparse(I, J, rand(rng_1, M)) # save in sparse format

# a_i Δcost/Δdemand
a_i = sparse(I, J, rand(rng_2, M)) # indexed by edges(SimpleDiGraph(A))

end # module
