"""
    TimeSeries <: TS

Composite type to store time series.
"""
struct TimeSeries{P <: SupportedTimeTypes} <: TS{P <: SupportedTimeTypes}
    times::AbstractVector{<:P}

    function TimeSeries(times)
        times = reshape(times, (length(times),))
        return new(times)
    end
end


"""
    ground_process(ts::TS)::TimeSeries

Returns the TimeSeries associated to a TS.
"""
function ground_process(ts::TS)::TimeSeries
    ts = TimeSeries(ts.times)
    return ts
end