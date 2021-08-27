module DynaPlanning

using LightGraphs, LightGraphsFlows
using GraphPlot, Plots
using Colors
using Compose, Cairo
using SparseArrays
using DataStructures
using Random
using Setfield
using LinearAlgebra
using JuMP, GLPK
using Logging

abstract type AbstractStrategy end

abstract type AbstractInstance end

abstract type AbstractAgent end

################################################################################

mutable struct EnvironmentModel
    model::Array{
        Pair{
            Tuple{LightGraphs.SimpleGraphs.SimpleEdge{Int64},AbstractStrategy},
            Tuple{Int64,LightGraphs.SimpleGraphs.SimpleEdge{Int64}},
        },
        1,
    }
    EnvironmentModel() = (x = new())
end

mutable struct StackelbergStrategy <: AbstractStrategy
    strategy::Array{Float64,1}
    StackelbergStrategy() = (x = new())
end

mutable struct SelfishStrategy <: AbstractStrategy
    strategy::Array{Float64,1}
    SelfishStrategy() = (x = new())
end

mutable struct Incentiver
    strategy_followed::AbstractStrategy
    Incentiver() = (x = new())
    Incentiver(AbstractStrategy) = new(AbstractStrategy)
end

mutable struct Q
    q::Array{
        Pair{Tuple{LightGraphs.SimpleGraphs.SimpleEdge{Int64},AbstractStrategy},Float64},
        1,
    }
    Q() = (x = new())
end

mutable struct DynaAgentPlus <: AbstractAgent
    path::Queue{Edge}
    state::Edge
    actions::Array{Edge,1}              # history
    γ::Float64                          # discount rate
    α::Float64                          # step-size
    ϵ::Float64                          # exploration parameter
    episodes::Int64                     # episodes
    τ::Int64                            # revisit time steps (Dyna-Q+)
    κ::Float64                          # weight
    θ::Float64                          # priority sweeping
    queue::PriorityQueue{Tuple{LightGraphs.SimpleGraphs.SimpleEdge{Int64},AbstractStrategy}}
    predecessors::Array{
        Pair{
            LightGraphs.SimpleGraphs.SimpleEdge{Int64},
            Array{Tuple{Int64,LightGraphs.SimpleGraphs.SimpleEdge{Int64},AbstractStrategy}},
        },
    }
    model::EnvironmentModel
    q::Q
    strategy::AbstractStrategy
    time::Int64
    DynaAgentPlus() = (x = new())
end

mutable struct StackelbergInstance <: AbstractInstance
    network::SimpleDiGraph
    costfunct::Array{Function,1}
    β::Float32                              # fraction
    agents::Array{DynaAgentPlus,1}
    StackelbergInstance() = (x = new())
end

end # module
