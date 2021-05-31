@testset "TS.jl" begin
    t1 = 1
    t2 = 2
    times = [t1, t2]
    marks = [4.2, -7]
    ts = SEMPP.TimeSeries(times)
    mts = SEMPP.MarkedTimeSeries(times, marks)
    times2 = [1, 2, 3]
    @test_throws ErrorException SEMPP.MarkedTimeSeries(times2, marks)
    @test SEMPP.ground_process(mts).times == ts.times
    @test SEMPP.start_time(ts) == t1
    @test SEMPP.start_time(mts) == t1
    @test SEMPP.end_time(ts) == t2
    @test SEMPP.end_time(mts) == t2
end