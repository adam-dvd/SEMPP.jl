using SEMPP
using Test, Random, Distributions, EGPD, CSV, DataFrames
using Dates

Random.seed!(12)

@testset "SEMPP.jl" begin
    include("Structures_tests.jl")
    include("ParametersEstimation_tests.jl")
    include(joinpath("Plots_tests", "PlotsData_tests.jl"))
end
