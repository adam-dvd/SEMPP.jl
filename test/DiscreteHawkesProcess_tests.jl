@testset "DiscreteHawkesProcess.jl" begin
    when = [1,2,3]
    times = [1,3]
    marks = [42.0, 4.2]
    pp = SEMPP.PointProcess(times)
    mpp = SEMPP.MarkedPointProcess(times, marks)
    μ = 3
    ϕ = 2.1
    γ = 7
    δ = 1.0
    @test_throws ErrorException SEMPP.volfunc(when, pp, -γ, δ)
    @test_logs (:warn,"no marks but δ non zero") SEMPP.volfunc(when, pp, γ, δ)
    @test size(SEMPP.volfunc(when, mpp, γ, δ)) == size(when)
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