ENV["JULIA_DEBUG"] = "DynaPlanning"

include("../src/DynaPlanning.jl")
include("../src/DynaUtils.jl")
include("TestDynaPlanning_Instances.jl")

using LightGraphs, LightGraphsFlows, GraphPlot, Compose # to reduce scope
using SparseArrays                      # idem
using DataStructures                    # idem
using Random
using Test
using Logging
using .Threads

import .DynaPlanning
import .DynaUtils
import .TestDynaPlanning_Instances

# Logging
logger = LogLevel(Logging.Debug)

################################################################################

@testset "Constructors" begin
    @test typeof(DynaPlanning.EnvironmentModel()) == Main.DynaPlanning.EnvironmentModel

    @test typeof(DynaPlanning.Incentiver()) == Main.DynaPlanning.Incentiver
    @test typeof(
        DynaPlanning.Incentiver(Main.DynaPlanning.StackelbergStrategy()).strategy_followed,
    ) == Main.DynaPlanning.StackelbergStrategy
    @test typeof(
        DynaPlanning.Incentiver(Main.DynaPlanning.SelfishStrategy()).strategy_followed,
    ) == Main.DynaPlanning.SelfishStrategy

    @test typeof(DynaPlanning.Q()) == Main.DynaPlanning.Q

    @test typeof(DynaPlanning.DynaAgentPlus()) == Main.DynaPlanning.DynaAgentPlus

    @test typeof(DynaPlanning.StackelbergInstance()) ==
          Main.DynaPlanning.StackelbergInstance
end

@testset "DynaUtils" begin
    ρ = DynaUtils.price_of_anarchy([0.0, 0.0, 1.0, 0.0], TestDynaPlanning_Instances.c_e)
    @test ρ == 4 / 3
end

@testset "Stackelberg" begin
    @test nthreads() == TestDynaPlanning_Instances.N_DynaAgents

    # find_flow in first epoch
    f_e_test = DynaUtils.find_flow(TestDynaPlanning_Instances.testStackelbergInstance)

    @test maximum(f_e_test) == 0.05
end
