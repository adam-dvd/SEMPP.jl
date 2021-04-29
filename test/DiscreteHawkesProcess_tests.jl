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

    @testset "negloglik(pp ; μ, ϕ, γ)" begin
        @test_logs (:warn, "ϕ must be positive or zero, taking absolute value") SEMPP.negloglik(pp, μ = μ, ϕ = -ϕ, γ = γ)
        @test SEMPP.negloglik(pp, μ = μ, ϕ = ϕ, γ = γ) isa Real
    end

    sepp = SEMPP.DiscreteSEPPExpKern(μ, ϕ, γ)
    
    @testset "negloglik(pp, sepp)" begin
        @test SEMPP.negloglik(pp, sepp) == SEMPP.negloglik(pp, μ = μ, ϕ = ϕ, γ = γ)
    end

    @testset "negloglik(mpp, mardens ; μ, ϕ, γ, δ, ξ, α, β, κ)" begin
        @test_logs (:warn, "μ γ must be positive or zero, taking absolute value") SEMPP.negloglik(mpp, GPD, μ = -μ, ϕ = ϕ, γ = -γ, ξ = ξ, α = α, β = β)
        @test SEMPP.negloglik(mpp, GPD, μ = -μ, ϕ = ϕ, γ = -γ, ξ = ξ, α = α, β = β) isa Real
    end

    sempp = SEMPP.DiscreteSEMPPExpKern(μ, ϕ, γ, δ, EGPD1, ξ, α, β, κ)

    @testset "negloglik(mpp, sempp)" begin
        @test SEMPP.negloglik(mpp, sempp) == SEMPP.negloglik(mpp, EGPD1, μ = μ, ϕ = ϕ, γ = γ, δ = δ, ξ = ξ, α = α, β = β, κ = κ)
    end

end