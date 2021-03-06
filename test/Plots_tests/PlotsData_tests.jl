@testset "PlotsData_tests.jl" begin
    when = [1,2,3]
    times = [1,3]
    marks = [2, 3]

    ts = SEMPP.TimeSeries(times)
    mts = SEMPP.MarkedTimeSeries(times, marks)

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
        s, p = SEMPP.pp_analysis(sepp, ts)
        @test length(s) == 2
        @test length(p) == 2

        s, p = SEMPP.pp_analysis(sempp_gpd, mts)
        @test length(s) == 2
        @test length(p) == 2
    end

    @testset "transformed_marks_ecdf" begin
        s, p = SEMPP.transformed_marks_ecdf(sempp_egpd, mts)
        @test length(s) == 2
        @test length(p) == 2

        s, p = SEMPP.transformed_marks_ecdf(sempp_gpd, mts)
        @test length(s) == 2
        @test length(p) == 2
    end

    @testset "marks_unit_exponential_qq" begin
        s, p = SEMPP.marks_unit_exponential_qq(sempp_egpd, mts)
        @test length(s) == 2
        @test length(p) == 2

        s, p = SEMPP.marks_unit_exponential_qq(sempp_gpd, mts)
        @test length(s) == 2
        @test length(p) == 2
    end
end