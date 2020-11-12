module DynaPlanning

using LightGraphs, LightGraphsFlows, SimpleWeightedGraphs
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

mutable struct LearningModel
    model::Array{
        Pair{
            Tuple{LightGraphs.SimpleGraphs.SimpleEdge{Int64},AbstractStrategy},
            Tuple{Int64,LightGraphs.SimpleGraphs.SimpleEdge{Int64}},
        },
        1,
    }
    LearningModel() = (x = new())
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
        Pair{
            Tuple{LightGraphs.SimpleGraphs.SimpleEdge{Int64},AbstractStrategy},
            Float64,
        },
        1,
    }
    Q() = (x = new())
end

mutable struct Dyna_agentPlus <: AbstractAgent
    path::Queue{Edge}
    state::Edge
    actions::Array{Edge,1}              # history
    γ::Float64                          # discount rate
    α::Float64                          # step-size
    ϵ::Float64                          # exploration parameter
    episodes::Int64
    τ::Int64                            # revisit time steps (Dyna-Q+)
    κ::Float64                          # weight
    θ::Float64                          # priority sweeping
    queue::PriorityQueue{
        Tuple{LightGraphs.SimpleGraphs.SimpleEdge{Int64},AbstractStrategy},
    }
    predecessors::Array{
        Pair{
            LightGraphs.SimpleGraphs.SimpleEdge{Int64},
            Array{
                Tuple{
                    Int64,
                    LightGraphs.SimpleGraphs.SimpleEdge{Int64},
                    AbstractStrategy,
                },
            },
        },
    }
    model::LearningModel
    q::Q
    strategy::AbstractStrategy
    time::Int64
    Dyna_agentPlus() = (x = new())
end

mutable struct StackelbergInstance <: AbstractInstance
    network::SimpleWeightedDiGraph
    β::Float32                              # fraction
    agents::Array{Dyna_agentPlus,1}
    StackelbergInstance() = (x = new())
end

################################################################################

"""
    find_flow(SInstance)

Compute the edge flow for `SInstance`.

"""
function find_flow(SInstance)

    W = LightGraphs.weights(SInstance.network)
    I = findnz(W.parent)[1]
    J = findnz(W.parent)[2]

    # init f_e
    f_e = sparse(I, J, zeros(ne(SInstance.network)))

    for a in SInstance.agents
        thisAgent = a
        thisPath = thisAgent.path
        edge_k = first(thisPath)
        f_e[edge_k.dst, edge_k.src] += 1.0
    end

    # normalize
    f_e = sparse(I, J, findnz(f_e)[3] / sum(findnz(f_e)[3]))

    return f_e
end

"""
    find_supply_demand(SInstance)

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
    total_demand = sum( [ x ] for x in d if x > zero(x) )

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

"""
    min_cost_flow(SInstance)

Compute the minimum-cost flow for `SInstance` using LP.

"""
function min_cost_flow(SInstance)

    # edges and capacity
    edgs = collect(edges(SInstance.network))

    # to Dict()
    c_dict = Dict()
    u_dict = Dict()
    edges_tuple = Array{Tuple{Int64,Int64},1}()

    for i in collect(1:ne(SInstance.network))
      c_dict[(edgs[i].src, edgs[i].dst)] = edgs[i].weight
      u_dict[(edgs[i].src, edgs[i].dst)] = 1.0
      push!(edges_tuple, (edgs[i].src, edgs[i].dst))
    end

    # supply / demand
    b = find_supply_demand( SInstance )

    nodes = collect(1:nv( SInstance.network ))

    min_costf = Model(GLPK.Optimizer)

    @variable(min_costf, 0 <= x[e in edges_tuple] <= u_dict[e])

    @objective(min_costf, Min, sum( c_dict[e] * x[e] for e in edges_tuple)  )

    for i in nodes
      @constraint(min_costf, sum(x[(ii,j)] for (ii,j) in edges_tuple if ii==i )
                      - sum(x[(j,ii)] for (j,ii) in edges_tuple if ii==i ) == b[i])
    end

    status = optimize!(min_costf)

    obj = objective_value(min_costf)

    x_star = value.(x)

    return x_star

end

function compute_Stackelberg_strategy(SInstance)

    x_star = min_cost_flow(SInstance)

    W = LightGraphs.weights(SInstance.network)

    I = findnz(W.parent)[1]
    J = findnz(W.parent)[2]

    m = sortperm(findnz(W.parent)[3])

    l = sparse( I[m], J[m], findnz(W.parent)[3][m] )

    βr = SInstance.β

    return strategies
end

function strategy_reward(incentiver, dyna_agent)
    if (isa(typeof(incentiver.strategy_followed), StackelbergStrategy))
        return 1
    else
        return -1
    end
end

function take_action(dyna_agent, G, Q)
    _reward = -999
    src = dyna_agent.state.src
    nghbrs = collect(neighbors(G, src))
    A = [LightGraphs.Edge(src, nghbrs[i]) for i = 1:length(nghbrs)]             # A
    S = dyna_agent.state                                                        # S ← current (nonterminal) state
    if rand(Float64, 1)[1] <= dyna_agent.ϵ                                      # A ← ϵ-greedy(S,Q)
        action = A[rand(1:length(A))]
    else
        for a = 1:length(A)
            idx_next_reward =
                findall(x -> x == (S, A[a]), broadcast(first, Q.q))
            if (!isempty(idx_next_reward))
                next_reward = Q.q[idx_next_reward].second
                if next_reward >= _reward
                    _reward = next_reward
                    action = A[rand(1:length(A))] # FIX ME
                end
            else
                action = A[rand(1:length(A))]
            end
        end
    end
    return action
end

function update_model(dyna_agent, action, reward)
    idx_model = findall(
        x -> x == (dyna_agent.state, action),
        broadcast(first, dyna_agent.model),
    )
    dyna_agent.model[idx_model] = (reward, action)
end

function update_state(dyna_agent, next_state)
    dyna_agent.state = next_state
end

function update_time(dyna_agent)
    dyna_agent.time += 1
end

function get_model(dyna_agent, S, A)
    idx_model = collect(findall(
        x -> x == (S, A),
        broadcast(first, dyna_agent.model.m),
    ))[1]
    return dyna_agent.model.m[idx_model].second[1],
    dyna_agent.model.m[idx_model].second[2],
    dyna_agent.time
end

function learn(dyna_agent, G, QQ, incentiver)
    print("Learning...")
    for ep = 1:dyna_agent.episodes
        while (dyna_agent.state != dyna_agent.path[end])
            a = take_action(dyna_agent, G, QQ)                       # Take action A; observe resultant reward R and state S'
            push!(dyna_agent.actions, a)

            # FIX ME
            a = LightGraphs.Edge(5, 1)
            b = LightGraphs.Edge(1, 10)
            #

            reward = strategy_reward(incentiver, dyna_agent)

            # priority sweeping
            _theta =
                reward + maximum(QQ.q[idx_QQ_next_state].second) -
                QQ.q[idx_QQ].second

            if _theta > dyna_agent.theta
                dyna_agent.queue.pop!((-_theta) => (dyna_agent.state, a))
            end

            push!(
                dyna_agent.model.m,
                (dyna_agent.state, dyna_agent.strategy) => (reward, a),
            )

            pop!(dyna_agent.predecessors, (reward, a))

            # Q(S,A) ← Q(S,A) + α [ R + γ max_a(Q(S',a)) - Q(S,A)]
            idx_QQ = collect(findall(
                x -> x == (dyna_agent.state, a),
                broadcast(first, QQ.q),
            ))[1]
            idx_QQ_next_state =
                collect(findall(x -> x == (a, b), broadcast(first, QQ.q)))[1]
            updated_QQ =
                QQ.q[idx_QQ].second +
                dyna_agent.α * (
                    reward + maximum(QQ.q[idx_QQ_next_state].second) -
                    QQ.q[idx_QQ].second
                )

            @set QQ.q[idx_QQ].second = updated_QQ

            idx_model = collect(findall(
                x -> x == (dyna_agent.state, a),
                broadcast(first, dyna_agent.model.m),
            ))[1]
            @show idx_model
            @set dyna_agent.model.m[idx_model].second = (reward, b) # Model(S,A) ← R,S'
            update_state(dyna_agent, b)
            update_time(dyna_agent)

            # Loop repeat n times
            for n = 1:dyna_agent.episodes
                # S ← random previously observed state
                #idx_model = rand(1:length(dyna_agent.model.m))
                # FIX ME
                idx_model = 1
                S = dyna_agent.model.m[idx_model].first[1]

                # A ← random action previously taken in S
                #idx_A = rand(1:length(dyna_agent.model.m))
                # FIX ME
                idx_A = 1
                A = dyna_agent.model.m[idx_model].first[2]

                # R,S' ← Model(S,A)
                _reward, next_state, _time = get_model(dyna_agent, S, A)
                _reward += dyna_agent.κ * sqrt(dyna_agent.time - _time)

                # Q(S,A) ← Q(S,A) + α [ R + γ max_a(Q(S',a)) - Q(S,A)]
                idx_QQ =
                    collect(findall(x -> x == (S, A), broadcast(first, QQ.q)))[1]
                idx_QQ_next_state =
                    collect(findall(x -> x == (A, b), broadcast(first, QQ.q)))[1]
                updated_QQ =
                    QQ.q[idx_QQ].second +
                    dyna_agent.α * (
                        reward + maximum(QQ.q[idx_QQ_next_state].second) -
                        QQ.q[idx_QQ].second
                    )

                @set QQ.q[idx_QQ].second = updated_QQ
            end
        end
    end
end

function _util_plot_graph(G, locs_x, locs_y)
    A = adjacency_matrix(G)
    nodelabel = collect(1:nv(SimpleDiGraph(A)))
    draw(
        PNG("graph.png", 10cm, 10cm),
        gplot(SimpleDiGraph(A), locs_x, locs_y, nodelabel = nodelabel),
    )
end

end # module
