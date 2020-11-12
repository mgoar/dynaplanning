ENV["JULIA_DEBUG"]="DynaPlanning"

include("DynaPlanning.jl")
include("testInstances.jl")

using LightGraphs, LightGraphsFlows, SimpleWeightedGraphs, GraphPlot, Compose # to reduce scope
using SparseArrays                      # idem
using DataStructures                    # idem
using Random
using Test
using Logging

import .DynaPlanning

# Logging
logger = LogLevel(Logging.Debug)

################################################################################

@testset "Constructors" begin
    @test typeof(DynaPlanning.LearningModel()) == Main.DynaPlanning.LearningModel

    @test typeof(DynaPlanning.Incentiver()) == Main.DynaPlanning.Incentiver
    @test typeof(DynaPlanning.Incentiver(Main.DynaPlanning.StackelbergStrategy()).strategy_followed) == Main.DynaPlanning.StackelbergStrategy
    @test typeof(DynaPlanning.Incentiver(Main.DynaPlanning.SelfishStrategy()).strategy_followed) == Main.DynaPlanning.SelfishStrategy

    @test typeof(DynaPlanning.Q()) == Main.DynaPlanning.Q

    @test typeof(DynaPlanning.Dyna_agentPlus()) == Main.DynaPlanning.Dyna_agentPlus

    @test typeof(DynaPlanning.StackelbergInstance()) == Main.DynaPlanning.StackelbergInstance
end;

@testset "Stackelberg" begin
    W = LightGraphs.weights(Main.testInstances.testStackelbergInstance.network)
    I = findnz(W.parent)[1]
    J = findnz(W.parent)[2]
    f_e_test = DynaPlanning.find_flow(Main.testInstances.testStackelbergInstance)
    @test findnz(f_e_test)[3][2] == 0.0

    allTestPaths = DynaPlanning.find_all_paths(Main.testInstances.testStackelbergInstance.network, 3, 10)
    @test length(allTestPaths) == 2

    supply_demand = DynaPlanning.find_supply_demand(Main.testInstances.testStackelbergInstance)
    @test length(supply_demand) == LightGraphs.nv(Main.testInstances.testStackelbergInstance.network)
    @test supply_demand[3] == 0.5
    @test supply_demand[end] == -0.5

    x_star = DynaPlanning.min_cost_flow(Main.testInstances.testStackelbergInstance)
end;
