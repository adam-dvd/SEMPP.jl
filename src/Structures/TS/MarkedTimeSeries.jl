"""
    MarkedTimeSeries <: TS
    
A composite type to store time data associated with marks.
"""
struct MarkedTimeSeries <: TS
    times::AbstractVector
    marks::AbstractVector

    function MarkedTimeSeries(times, marks) 
        length(times) == length(marks) || error("times and marks must be the same size")
        times = reshape(times, (length(times),))
        marks = reshape(marks, (length(marks),))
        return new(times, marks)
    end
end

"""
    MarkedTimeSeries(ts::TS)::MarkedTimeSeries

A method to see a point process as marked with zero marks.
"""
function MarkedTimeSeries(ts::TS)::MarkedTimeSeries
    if ts isa TimeSeries
        times = ts.times
        marks = fill(0, size(times))
        return MarkedTimeSeries(times, marks)
    end
    
    return ts
end