using SEMPP
using Test, Random

Random.seed!(12)

@testset "SEMPP.jl" begin
    include("PP_tests.jl")
    include("DiscreteHawkesProcess_tests.jl")
end
