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

    sepp = SEMPP.SEPPExpKern(μ, ϕ, γ)
    
    @testset "discrete_negloglik(pp, sepp)" begin
        @test SEMPP.discrete_negloglik(pp, sepp) == SEMPP.discrete_negloglik(pp, μ = μ, ϕ = ϕ, γ = γ)
    end

    @testset "discrete_negloglik(mpp, mardens ; μ, ϕ, γ, δ, ξ, α, β, κ)" begin
        @test_logs (:warn, "μ γ must be positive or zero, taking absolute value") SEMPP.discrete_negloglik(mpp, GPD, μ = -μ, ϕ = ϕ, γ = -γ, ξ = ξ, α = α, β = β)
        @test SEMPP.discrete_negloglik(mpp, GPD, μ = -μ, ϕ = ϕ, γ = -γ, ξ = ξ, α = α, β = β) isa Real
    end

    sempp_egpd = SEMPPExpKern(μ, ϕ, γ, δ, EGPD1, ξ, α, β, κ)
    sempp_gpd = SEMPPExpKern(μ, ϕ, γ, δ, GPD, ξ, α, β)

    @testset "discrete_negloglik(mpp, sempp)" begin
        @test SEMPP.discrete_negloglik(mpp, sempp_egpd) == SEMPP.discrete_negloglik(mpp, EGPD1, μ = μ, ϕ = ϕ, γ = γ, δ = δ, ξ = ξ, α = α, β = β, κ = κ)
    end

    @testset "discrete_fit!(sepp, pp)" begin
        @test discrete_fit!(sepp, pp) isa Real
        @test μ != sepp.μ
        @test ϕ != sepp.ϕ
        @test γ != sepp.γ
        @test sepp.μ >= 0
        @test sepp.ϕ >= 0
        @test sepp.γ >= 0
    end

    @testset "discrete_fit!(sempp, mpp)" begin
        @test discrete_fit!(sempp_egpd, mpp) isa Real
        @test discrete_fit!(sempp_gpd, mpp) isa Real
        @test μ != sempp_egpd.μ
        @test ϕ != sempp_egpd.ϕ
        @test γ != sempp_egpd.γ
        @test δ != sempp_egpd.δ
        @test ξ != sempp_egpd.ξ
        @test α != sempp_egpd.α
        @test β != sempp_egpd.β
        @test κ != sempp_egpd.κ

        @test sempp_egpd.μ >=0
        @test sempp_egpd.ϕ >=0
        @test sempp_egpd.γ >=0
        @test sempp_egpd.δ >=0
        @test sempp_egpd.ξ >=0
        @test sempp_egpd.α >=0
        @test sempp_egpd.β >=0
        @test sempp_egpd.κ >=0

        @test μ != sempp_gpd.μ
        @test ϕ != sempp_gpd.ϕ
        @test γ != sempp_gpd.γ
        @test δ != sempp_gpd.δ
        @test ξ != sempp_gpd.ξ
        @test α != sempp_gpd.α
        @test β != sempp_gpd.β

        @test sempp_gpd.μ >=0
        @test sempp_gpd.ϕ >=0
        @test sempp_gpd.γ >=0
        @test sempp_gpd.δ >=0
        @test sempp_gpd.ξ >=0
        @test sempp_gpd.α >=0
        @test sempp_gpd.β >=0
    end

end