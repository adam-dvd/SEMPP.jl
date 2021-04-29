@testset "DiscreteHawkesProcess.jl" begin

    when = [1,2,3]
    times = [1,3]
    marks = [2, 3]

    pp = SEMPP.PointProcess(times)
    mpp = SEMPP.MarkedPointProcess(times, marks)

    μ = 0
    ϕ = 1
    γ = 1
    δ = 1.0

    @testset "volfunc" begin
        @test_throws ErrorException SEMPP.volfunc(when, pp, -γ, δ)
        @test_logs (:warn,"no marks but δ non zero") SEMPP.volfunc(when, pp, γ, δ)
        @test size(SEMPP.volfunc(when, mpp, γ, δ)) == size(when)
        @test SEMPP.volfunc(when, mpp, γ, δ) == [0, 3*exp(-1), 3*exp(-2)]
    end

    GPD = Distributions.GeneralizedPareto
    EGPD1 = EGPD.EGPpower
    ξ = 1.5
    β = 2
    α = 1.1
    κ = 1.5
    #=
    @test_warn negloglik(times, marks, GPD, μ = μ, ϕ = -ϕ, γ = γ, δ = δ, ξ = ξ, β = β, α = α, κ = κ)
    @test negloglik(times, marks, GPD, μ = μ, ϕ = ϕ, γ = γ, δ = δ, ξ = ξ, β = β, α = α, κ = κ) isa Real
    @test negloglik(times, marks, EGPD1, μ = μ, ϕ = ϕ, γ = γ, δ = δ, ξ = ξ, β = β, α = α, κ = κ) isa Real
    =#
end