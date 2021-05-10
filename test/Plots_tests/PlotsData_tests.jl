@testset "PlotsData_tests.jl" begin
    when = [1,2,3]
    times = [1,3]
    marks = [2, 3]

    pp = SEMPP.PointProcess(times)
    mpp = SEMPP.MarkedPointProcess(times, marks)

    μ = 2.0
    ϕ = 1
    γ = 1
    δ = 1.0

    GPD = Distributions.GeneralizedPareto
    EGPD1 = EGPD.EGPpower

    ξ = 1.5
    β = 2
    α = 1.1
    κ = 1.5

    sepp = SEMPP.SEPPExpKern(μ, ϕ, γ)
    sempp_egpd = SEMPPExpKern(μ, ϕ, γ, δ, EGPD1, ξ, α, β, κ)
    sempp_gpd = SEMPPExpKern(μ, ϕ, γ, δ, GPD, ξ, α, β)

    @testset "pp_analysis" begin
        s, p = pp_analysis(sepp, pp)
        @test length(s) == 2
        @test length(p) == 2

        s, p = pp_analysis(sempp, mpp)
        @test length(s_m) == 2
        @test length(p_m) == 2
    end

    @testset "transformed_marks_ecdf" begin
        s, p = transformed_marks_ecdf(sempp_egpd, mpp)
        @test length(s) == 2
        @test length(p) == 2

        s, p = transformed_marks_ecdf(sempp_gpd, mpp)
        @test length(s) == 2
        @test length(p) == 2
    end

    @testset "marks_unit_exponential_qq" begin
        s, p = marks_unit_exponential_qq(sempp_egpd, mpp)
        @test length(s) == 2
        @test length(p) == 2

        s, p = marks_unit_exponential_qq(sempp_gpd, mpp)
        @test length(s) == 2
        @test length(p) == 2
    end
end