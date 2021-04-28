using SEMPP

@testset "PP.jl" begin
    t1 = 1
    t2 = 2
    times = [t1, t2]
    marks = [4.2, -7]
    pp = SEMPP.PointProcess(times)
    mpp = SEMPP.MarkedPointProcess(times, marks)
    @test SEMPP.ground_process(mpp).times == pp.times
    @test SEMPP.start_time(pp) == t1
    @test SEMPP.start_time(mpp) == t1
    @test SEMPP.end_time(pp) == t2
    @test SEMPP.end_time(mpp) == t2
end