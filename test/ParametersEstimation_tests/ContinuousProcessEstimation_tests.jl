@testset "ProcessEstimation.jl" begin
    when = [1,2,3]
    times = [1,3]
    marks = [2, 3]

    ts = TimeSeries(times)
    mts = MarkedTimeSeries(times, marks)

    μ = 2.0
    ϕ = 1
    γ = 1
    δ = 1.0

    ξ = 1.5
    β = 2
    α = 1.1
    κ = 1.5

    @testset "negloglik(ts ; μ, ϕ, γ)" begin
        @test_logs (:warn, "ϕ must be positive or zero, taking absolute value") SEMPP.negloglik(ts, μ = μ, ϕ = -ϕ, γ = γ)
        @test SEMPP.negloglik(ts, μ = μ, ϕ = ϕ, γ = γ) isa Real
    end

    sepp = SEMPP.SEPPExpKern(ts, μ = μ, ϕ = ϕ, γ = γ)
    
    @testset "negloglik(sepp)" begin
        @test SEMPP.negloglik(sepp) == SEMPP.negloglik(ts, μ = μ, ϕ = ϕ, γ = γ)
    end

    @testset "negloglik(mts, mardens ; μ, ϕ, γ, δ, ξ, α, β, κ)" begin
        @test_logs (:warn, "μ γ must be positive or zero, taking absolute value") SEMPP.negloglik(mts, Distributions.GeneralizedPareto, μ = -μ, ϕ = ϕ, γ = -γ, ξ = ξ, α = α, β = β)
        @test SEMPP.negloglik(mts, Distributions.GeneralizedPareto, μ = -μ, ϕ = ϕ, γ = -γ, ξ = ξ, α = α, β = β) isa Real
    end

    sempp_egpd = SEMPPExpKern(mts, μ = μ, ϕ = ϕ, γ = γ, δ = δ, markdens = EGPD.EGPpower, ξ = ξ, α = α, β = β, κ = κ)
    sempp_gpd = SEMPPExpKern(mts, μ = μ, ϕ = ϕ, γ = γ, δ = δ, markdens = Distributions.GeneralizedPareto, ξ = ξ, α = α, β = β)

    @testset "negloglik(sempp)" begin
        @test SEMPP.negloglik(sempp_egpd) == SEMPP.negloglik(mts, EGPD.EGPpower, μ = μ, ϕ = ϕ, γ = γ, δ = δ, ξ = ξ, α = α, β = β, κ = κ)
    end

    @testset "fit!(sepp)" begin
        @test fit!(sepp) isa Real
        @test μ != sepp.μ
        @test ϕ != sepp.ϕ
        @test γ != sepp.γ
        @test sepp.μ >= 0
        @test sepp.ϕ >= 0
        @test sepp.γ >= 0
    end

    @testset "fit!(sempp)" begin
        @test fit!(sempp_egpd) isa Real
        @test fit!(sempp_gpd) isa Real
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

    @testset "Simulation.jl" begin

        @testset "PointProcess" begin
            simulated_ts = simulation(sepp, 100*365)
        
            @test simulated_ts isa TimeSeries

            sepp_simul = SEPPExpKern(simulated_ts)
            fit!(sepp_simul)

            @test abs(sepp_simul.μ - sepp.μ) < 1.96 * sepp_simul.cov_mat[1,1]
            @test abs(sepp_simul.ϕ - sepp.ϕ) < 1.96 * sepp_simul.cov_mat[2,2]
            @test abs(sepp_simul.γ - sepp.γ) < 1.96 * sepp_simul.cov_mat[3,3]
        end
        
        @testset "MarkedPointProcess GPD" begin
            simulated_mts = simulation(sempp_gpd, 100*365)

            @test simulated_mts isa MarkedTimeSeries

            sempp_gpd_simul = SEMPPExpKern(simulated_mts)
            fit!(sempp_gpd_simul)

            @test abs(sempp_gpd_simul.μ - sempp_gpd.μ) < 1.96 * sempp_gpd_simul.cov_mat[1,1]
            @test abs(sempp_gpd_simul.ϕ - smpp_gpd.ϕ) < 1.96 * sempp_gpd_simul.cov_mat[2,2]
            @test abs(sempp_gpd_simul.γ - sempp_gpd.γ) < 1.96 * sempp_gpd_simul.cov_mat[3,3]
        end
        
        @testset "MarkedPointProcess EGPD" begin
            simulated_mts = simulation(sempp_egpd, 100*365)

            @test simulated_mts isa MarkedTimeSeries

            sempp_egpd_simul = SEMPPExpKern(simulated_mts)
            fit!(sempp_egpd_simul)

            @test abs(sempp_egpd_simul.μ - sempp_egpd.μ) < 1.96 * sempp_egpd_simul.cov_mat[1,1]
            @test abs(sempp_egpd_simul.ϕ - smpp_egpd.ϕ) < 1.96 * sempp_egpd_simul.cov_mat[2,2]
            @test abs(sempp_egpd_simul.γ - sempp_egpd.γ) < 1.96 * sempp_egpd_simul.cov_mat[3,3]    
        end
    end
end