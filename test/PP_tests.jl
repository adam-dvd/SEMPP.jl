@testset "PP.jl" begin
    t1 = 1
    t2 = 2
    times = [t1, t2]
    marks = [4.2, -7]
    pp = SEMPP.PointProcess(times)
    mpp = SEMPP.MarkedPointProcess(times, marks)
    times2 = [1, 2, 3]
    @test_throws ErrorException SEMPP.MarkedPointProcess(times2, marks)
    @test SEMPP.ground_process(mpp).times == pp.times
    @test SEMPP.start_time(pp) == t1
    @test SEMPP.start_time(mpp) == t1
    @test SEMPP.end_time(pp) == t2
    @test SEMPP.end_time(mpp) == t2
end