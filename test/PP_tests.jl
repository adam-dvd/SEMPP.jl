@testset "PP.jl" begin
    t1 = DateTime(2021,04,21,18,36)
    t2 = DateTime(2021,04,21,18,37)
    times = [t1, t2]
    marks = [4.2, -7]
    pp = PointProcess(times)
    mpp = MarkedPointProcess(times, marks)
    @test ground_process(mpp).times == pp.times
    @test start_time(pp) == t1
    @test start_time(mpp) == t1
    @test end_time(pp) == t2
    @test end_time(mpp) == t2
end