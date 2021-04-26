"""
A type for point process might be useful and precise things, but it might conversely impede clarity... not sure to keep it
"""

abstract type PP end


struct PointProcess <: PP
    times::AbstractVector
end


struct MarkedPointProcess <: PP
    times::AbstractVector
    marks::AbstractVector
    MarkedPointProcess(times, marks) = size(times) == size(marks) ? new(times, marks) : error("times and marks must be the same size")
end


function ground_process(mpp::MarkedPointProcess)
    pp = PointProcess(mpp.times)
    return pp
end


function start_time(pp::PP)
    return first(pp.times)
end


function end_time(pp::PP)
    return last(pp.times)
end

