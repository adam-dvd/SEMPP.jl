"""
    volfunc(when::AbstractVector, ts::TS, γ::Real, δ::Real = 0)::AbstractVector{<:Real}

Compute the volatility function that is the term in the rate of a SEPP that corresponds to the self-excitement, denoted by ν in Li2020.
"""
function volfunc(when::AbstractVector, ts::TS, γ::Real, δ::Real = 0)::AbstractVector{<:Real}
    
    γ < 0 && error("γ must be positive or zero")

    (ts isa TimeSeries && δ != 0) && (@warn "no marks but δ non zero")

    first(when) isa TimeType && (when = Dates.value.(when))
    times = first(ts.times) isa TimeType ? Dates.value.(ts.times) : ts.times
    marks = ts isa MarkedTimeSeries ? ts.marks : fill(0, size(times))

    mts = hcat(times, marks)


    function self_ex(t)

        mts_to_t = mts[times .< t, :]

        term(markedpoint) = (1+ δ*markedpoint[2])*exp(-γ*(t-markedpoint[1]))

        return sum(term.(eachrow(mts_to_t)))
    end


    return self_ex.(when)
end

include("ParametersEstimation/DiscreteProcessEstimation.jl")
include("ParametersEstimation/ContinuousProcessEstimation.jl")