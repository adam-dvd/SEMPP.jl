@testset "ParametersEstimation.jl" begin
    when = [1,2,3]
    times = [1,3]
    marks = [2, 3]

    ts = TimeSeries(times)
    mts = MarkedTimeSeries(times, marks)

    μ = 2.0
    ϕ = 1
    γ = 1
    δ = 1.0

    @testset "volfunc" begin
        @test_throws ErrorException SEMPP.volfunc(when, ts, -γ, δ)
        @test_logs (:warn,"no marks but δ non zero") SEMPP.volfunc(when, ts, γ, δ)
        @test size(SEMPP.volfunc(when, mts, γ, δ)) == size(when)
        @test SEMPP.volfunc(when, mts, γ, δ) == [0, 3*exp(-1), 3*exp(-2)]
    end
end

include("ParametersEstimation_tests/DiscreteProcessEstimation_tests.jl")
include("ParametersEstimation_tests/ContinuousProcessEstimation_tests.jl")