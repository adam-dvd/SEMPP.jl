using SEMPP
using Test, Random, Distributions, EGPD

Random.seed!(12)

@testset "SEMPP.jl" begin
    include("Types_tests.jl")
    include("DiscreteHawkesProcess_tests.jl")
end
