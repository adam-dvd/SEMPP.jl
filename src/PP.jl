"""
A type for point process might be useful and precise things, but it might conversely impede clarity... not sure to keep it
"""

abstract type PP end


struct PointProcess <: PP
    times::Array{TimeType}
end


struct MarkedPointProcess <: PP
    times::Array{TimeType}
    marks::Array{Real}
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

