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
using .Threads
using Actors
import Actors: spawn

abstract type AbstractStrategy end

abstract type AbstractInstance end

abstract type AbstractAgent end

################################################################################

mutable struct DynaTuple
    e::LightGraphs.SimpleGraphs.SimpleEdge{Int64}
    strategy::AbstractStrategy
    s::LightGraphs.SimpleGraphs.SimpleEdge{Int64}
end

mutable struct DictSrv{L}
    lk::L
end

(ds::DictSrv)(f::Function, args...) = call(ds.lk, f, args...)

# Indexing interface
Base.getindex(d::DictSrv, key) = call(d.lk, getindex, key)
Base.setindex!(d::DictSrv, value, key) = call(d.lk, setindex!, value, key)

# Behavior
ds(d::Dict, f::Function, args...) = f(d, args...)

# Start dictsrv
dictsrv(d::Dict; remote = false) = DictSrv(spawn(ds, d; remote))

export DictSrv, dictsrv

mutable struct EnvironmentModel
    model::DictSrv
    EnvironmentModel() = (x = new(); x.model = dictsrv(Dict{DynaTuple,Int}()))

    function Base.isequal(A::DynaTuple, B::DynaTuple)
        A.e == B.e && A.strategy == B.strategy && A.s == B.s
    end

    function Base.hash(A::DynaTuple)
        hash(A.e + A.strategy + A.s)
    end
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
        Pair{LightGraphs.SimpleGraphs.SimpleEdge{Int64},AbstractStrategy},
        Pair{LightGraphs.SimpleGraphs.SimpleEdge{Int64},Float64},
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

    function chooseAction(a)

        max_reward = -9999
        # ϵ-greedy
        if rand(Float64, 1) <= ϵ
            # randomly choose action
            action = a[rand(1:end)]
        else
            # greedy action
            current_state = state
            if length(Set(q)) == 1
                action = rand(actions)
            else
                for a in actions
                    this_reward = q[current_state][a]
                    if this_reward >= max_reward
                        action = a
                        max_reward = this_reward
                    end
                end
            end
        end

        return action

    end

    function learn()



    end
end

mutable struct StackelbergInstance <: AbstractInstance
    network::SimpleDiGraph
    costfunct::Array{Function,1}
    β::Float32                              # fraction
    agents::Array{DynaAgentPlus,1}
    StackelbergInstance() = (x = new())
end

end # module
