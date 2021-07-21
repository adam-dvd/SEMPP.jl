"""
    discrete_tail_estimation(sempp::SEMPPExpKern, q::Real = 0.95; to_time::Union{DiscreteTimeTypes, Nothing} = nothing, from_time::Union{DiscreteTimeTypes, Nothing} = nothing)::Real

Compute the quantile at level *q* of the marks distribution between from_time and to_time, knowing history until from_time, with a discrete time approach.

See Li2020 for the mathematics behind this method.

from_time defaults to the last timestamp in the time series, and to_time defaults to from_time + 1. 
Note that if from_time differs from the last timestamp in the time series, it will be considered that no event happened between this last timestamp and from_time.
"""
function discrete_tail_estimation(sempp::SEMPPExpKern, q::Real = 0.95; to_time::Union{DiscreteTimeTypes, Nothing} = nothing, from_time::Union{DiscreteTimeTypes, Nothing} = nothing, history::Union{MarkedTimeSeries, Nothing} = nothing)::Real 
    mts = isnothing(history) ? MarkedTimeSeries([], []) : history
    isnothing(mts) && error("No data in model, can't estimate")

    isnothing(from_time) && (from_time = end_time(mts))

    if start_time(mts) isa TimeType

        isnothing(to_time) && (to_time = from_time + Day(1))

        if !(from_time isa TimeType) || !(to_time isa TimeType) 
            error("times and from_date and to_date must be either all Integers or all TimeTypes")
        end

        when = from_time:Dates.Day(1):to_time
    else

        isnothing(to_time) && (to_time = from_time + 1)

        if !(from_time isa Integer) || !(to_time isa Integer) 
            error("times and from_date and to_date must be either all Integers or all TimeTypes")
        end

        when = from_time:to_time
    end

    μ = sempp.μ
    ϕ = sempp.ϕ
    γ = sempp.γ
    δ = sempp.δ

    vol = volfunc(when, mts, γ, δ)
    intens = sum(μ .+ ϕ .* vol)
    prob = 1 - exp(-intens )

    prob <= 1-q && (error("Probability of an event below $(1-q), cannot estimate tail."))

    markdens = sempp.markdens
    ξ = sempp.ξ
    α = sempp.α
    β = sempp.β
    κ = sempp.κ

    σ = β + α * vol

    if markdens == Distributions.GeneralizedPareto
        exceed_q = σ / ξ * ((prob / (1 - q))^ξ - 1)
    else
        exceed_q = σ / ξ * ((1 - (1 - (1 - q) / prob)^(1/κ))^(-ξ) - 1)
    end

    return exceed_q
end


"""
    tail_estimation(sempp::SEMPPExpKern, q::Real = 0.95; to_time::Union{DiscreteTimeTypes, Nothing} = nothing, from_time::Union{DiscreteTimeTypes, Nothing} = nothing)::Real

Compute the quantile at level *q* of the marks distribution between from_time and to_time, knowing history until from_time, with a continuous time approach.

See Li2020 for the mathematics behind this method.

from_time defaults to the last timestamp in the time series, and to_time defaults to from_time + 1. 
Note that if from_time differs from the last timestamp in the time series, it will be considered that no event happened between this last timestamp and from_time.
"""
function tail_estimation(sempp::SEMPPExpKern, q::Real = 0.95; to_time::Union{DiscreteTimeTypes, Nothing} = nothing, from_time::Union{DiscreteTimeTypes, Nothing} = nothing)::Real 
    mts = sempp.data
    isnothing(mts) && error("No data in model, can't estimate")

    times = mts.times
    marks = mts.marks

    isnothing(from_time) && (from_time = last(times))

    if first(times) isa TimeType

        isnothing(to_time) && (to_time = from_time + Day(1))

        if !(from_time isa TimeType) || !(to_time isa TimeType) 
            error("times and from_date and to_date must be either all Real or all TimeTypes")
        end

        Dates.value(Date(to_time - from_time))

    else

        isnothing(to_time) && (to_time = from_time + 1)

        if (from_time isa TimeType) || (to_time isa TimeType) 
            error("times and from_date and to_date must be either all Real or all TimeTypes")
        end

        T = to_time - from_time

    end

    μ = sempp.μ
    ϕ = sempp.ϕ
    γ = sempp.γ
    δ = sempp.δ

    intens = μ * T + ϕ/γ * sum((1 .+ δ .* marks) .* (1 .- exp.(-γ .* (to_time .- times))))
    prob = 1 - exp(-intens )

    prob <= 1-q && (error("Probability of an event below $(1-q), cannot estimate tail."))

    markdens = sempp.markdens
    ξ = sempp.ξ
    α = sempp.α
    β = sempp.β
    κ = sempp.κ

    σ = β + α * vol

    if markdens == Distributions.GeneralizedPareto
        exceed_q = σ / ξ * ((prob / (1 - q))^ξ - 1)
    else
        exceed_q = σ / ξ * ((1 - (1 - (1 - q) / prob)^(1/κ))^(-ξ) - 1)
    end

    return exceed_q
end