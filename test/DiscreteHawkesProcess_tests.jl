@testset "DiscreteHawkesProcess.jl" begin

    when = [1,2,3]
    times = [1,3]
    marks = [2, 3]

    pp = SEMPP.PointProcess(times)
    mpp = SEMPP.MarkedPointProcess(times, marks)

    μ = 2.0
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

    @testset "discrete_negloglik(pp ; μ, ϕ, γ)" begin
        @test_logs (:warn, "ϕ must be positive or zero, taking absolute value") SEMPP.discrete_negloglik(pp, μ = μ, ϕ = -ϕ, γ = γ)
        @test SEMPP.discrete_negloglik(pp, μ = μ, ϕ = ϕ, γ = γ) isa Real
    end

    sepp = SEMPP.DiscreteSEPPExpKern(μ, ϕ, γ)
    
    @testset "discrete_negloglik(pp, sepp)" begin
        @test SEMPP.discrete_negloglik(pp, sepp) == SEMPP.discrete_negloglik(pp, μ = μ, ϕ = ϕ, γ = γ)
    end

    @testset "discrete_negloglik(mpp, mardens ; μ, ϕ, γ, δ, ξ, α, β, κ)" begin
        @test_logs (:warn, "μ γ must be positive or zero, taking absolute value") SEMPP.discrete_negloglik(mpp, GPD, μ = -μ, ϕ = ϕ, γ = -γ, ξ = ξ, α = α, β = β)
        @test SEMPP.discrete_negloglik(mpp, GPD, μ = -μ, ϕ = ϕ, γ = -γ, ξ = ξ, α = α, β = β) isa Real
    end

    sempp = SEMPP.DiscreteSEMPPExpKern(μ, ϕ, γ, δ, EGPD1, ξ, α, β, κ)

    @testset "discrete_negloglik(mpp, sempp)" begin
        @test SEMPP.discrete_negloglik(mpp, sempp) == SEMPP.discrete_negloglik(mpp, EGPD1, μ = μ, ϕ = ϕ, γ = γ, δ = δ, ξ = ξ, α = α, β = β, κ = κ)
    end

    @testset "fit!(sepp, pp)" begin
        @test SEMPP.fit!(sepp, pp) isa Real
        @test μ != sepp.μ
        @test ϕ != sepp.ϕ
        @test γ != sepp.γ
        @test sepp.μ >= 0
        @test sepp.ϕ >= 0
        @test sepp.γ >= 0

    end

end