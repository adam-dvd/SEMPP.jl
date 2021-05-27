"""
    DiscreteTimeTypes

Type for discrete time handling.
"""
DiscreteTimeTypes = Union{Integer, TimeType}

"""
    SupportedTimeTypes

Types supported for the timestamps data in the time series.
"""
SupportedTimeTypes = Union{Real, DiscreteTimeTypes}


"""
    TS{P}

Abstract type to store data. P is the type of the time stamps in the stored time series.

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