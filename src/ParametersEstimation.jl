"""
    volfunc(when::AbstractVector, ts::TS, γ::Real, δ::Real = 0)::AbstractVector{<:Real}

Compute the volatility function that is the term in the rate of a SEPP that corresponds to the self-excitement, denoted by ν in Li2020.
"""
function volfunc(when::AbstractVector, ts::TS, γ::Real, δ::Real = 0, impact_func::Function = (d -> 1 + tanh(d)))::AbstractVector{<:Real}

    isempty(when) && (return when)
    isempty(ts.times) && return(zero(when))
    
    γ < 0 && error("γ must be positive or zero")

    (ts isa TimeSeries && δ != 0) && (@warn "no marks but δ non zero")

    first(when) isa TimeType && (when = Dates.value.(DateTime.(when)) ./ (1000*3600*24))
    times = first(ts.times) isa TimeType ? Dates.value.(DateTime.(ts.times)) ./ (1000*3600*24) : copy(ts.times)
    marks = ts isa MarkedTimeSeries ? copy(ts.marks) : fill(0, size(times))

    mts = hcat(times, marks)


    function self_ex(t)

        mts_to_t = mts[times .< t, :]

        term(markedpoint) = impact_func(δ * markedpoint[2]) * exp(-γ * (t - markedpoint[1]))

        isempty(mts_to_t) && return 0

        return sum(term.(eachrow(mts_to_t)))
    end


    return self_ex.(when)
end

include("ParametersEstimation/DiscreteProcessEstimation.jl")
include("ParametersEstimation/ContinuousProcessEstimation.jl")