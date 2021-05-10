using SEMPP
using Test, Random, Distributions, EGPD, DataFrames, Gadfly

Random.seed!(12)

@testset "SEMPP.jl" begin
    include("Types_tests.jl")
    include("DiscreteHawkesProcess_tests.jl")
    include(joinpath("Plots_tests", "PlotsData_tests.jl"))
end
