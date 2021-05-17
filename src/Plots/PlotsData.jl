"""
    ecdf(y::Vector{<:Real})::Tuple{Vector{<:Real}, Vector{<:Real}}

Compute the empirical cumulative distribution function using the Gumbel formula.

The empirical quantiles are computed using the Gumbel plotting positions as
as recommended by [Makkonen (2006)](https://journals.ametsoc.org/jamc/article/45/2/334/12668/Plotting-Positions-in-Extreme-Value-Analysis).

# Example
```julia-repl
julia> (x, F̂) = Extremes.ecdf(y)
```
# Reference
Makkonen, L. (2006). Plotting positions in extreme value analysis. Journal of
Applied Meteorology and Climatology, 45(2), 334-340.
"""
function ecdf(y::Vector{<:Real})::Tuple{Vector{<:Real}, Vector{<:Real}}
    ys = sort(y)
    n = length(ys)
    p = collect(1:n)/(n+1)

    return ys, p
end

"""
    pp_analysis(sepp::SEPP)

Compute the rescaled time transformed of the realizations of the point process.

It should behave like order statistics from a standard uniform distribution.
"""
function pp_analysis(sepp::SEPP)
    ts=sepp.data
    isnothing(ts) && error("No data in model")

    mts = MarkedTimeSeries(ts)
    marks = mts.marks
    times = first(mts.times) isa TimeType ? Dates.value.(mts.times) : mts.times

    starttime = first(times)

    μ = sepp.μ
    ϕ = sepp.ϕ
    γ = sepp.γ
    δ = sepp isa SEMPPExpKern ? sepp.δ : 0

    aux(i) = sum((1 .+ δ .* marks[1:i]) .* (1 .- exp.(-γ .*(times[i] .- times[1:i]))))

    τs = μ .* (times .- starttime) .+ ϕ/γ .* aux.(collect(1:length(times)))

    scaled = τs ./ (last(τs))

    n = length(scaled)

    return scaled, (1/n):(1/n):1
end

"""
    transformed_marks_ecdf(sempp::SEMPPExpKern)

Compute the transformed marks ecdf for ploting.

See Li2020 4.2..
"""
function transformed_marks_ecdf(sempp::SEMPPExpKern)
    mts = sempp.data
    isnothing(mts) && error("No data in model")
    
    times = mts.times
    marks = mts.marks
    
    γ = sempp.γ
    δ = sempp.δ
    markdens = sempp.markdens
    ξ = sempp.ξ
    α = sempp.α
    β = sempp.β
    κ = sempp.κ

    vol = volfunc(times, mts, γ, δ)
    σ = β .+ α .* vol

    sigmarks = hcat(σ, marks)

    if markdens == Distributions.GeneralizedPareto
        u = (sigmark -> cdf(Distributions.GeneralizedPareto(0, sigmark[1], ξ), sigmark[2])).(eachrow(sigmarks))
    else    # EGPD case
        u = (sigmark -> cdf(EGPD.EGPpower(sigmark[1], ξ, κ), sigmark[2])).(eachrow(sigmarks))
    end

    return ecdf(u)
end

"""
    marks_unit_exponential_qq(sempp::SEMPPExpKern)
Compute the quantile plot data based on the unit exponential distribution.

See 4.2 in Li2020.
"""
function marks_unit_exponential_qq(sempp::SEMPPExpKern)

    u, p = transformed_marks_ecdf(sempp)
    
    model_qs = -log.(1 .- p)

    emp_qs = -log.(1 .- u)

    return emp_qs, model_qs
end