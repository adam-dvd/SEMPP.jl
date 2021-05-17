"""
    TS

Abstract type to store data.

Subtypes are :
    - TimeSeries
    - MarkedTimeSeries
"""
abstract type TS end


"""
    start_time(ts::TS)

Get the start time of a time series.
"""
function start_time(ts::TS)
    return first(ts.times)
end


"""
    end_time(ts::TS)

Get the end time of a time series.
"""
function end_time(ts::TS)
    return last(ts.times)
end


include("TS/TimeSeries.jl")
include("TS/MarkedTimeSeries.jl")


