"""
    SupportedMarksDistributions

Union type for the marks distributions supported by the package for the SEMPP models.
"""
SupportedMarksDistributions = Union{Type{Distributions.GeneralizedPareto}, Type{EGPD.EGPpower}}


"""
SEMPPExpKern(data, μ, ϕ, γ, markdens, ξ, α, β, κ)

Self-exciting marked point process with exponential kernel model type such that for all t, and all m there is an event in t with mark m with probability λ_g(t)f(m|t)) where :
```
λ_g(t) = μ + ϕ * ν(t)
```

with `ν(t) = \\sum_{t_k < t} (1 + δ*m_k) (\\exp(-γ(t-t_k))`

f(m|t) is markdens with scale `σ_t = β + α * ν(t)`

`t_k` are the timestamps of all events and `m_k` their marks.
"""
mutable struct SEMPPExpKern <: SEPP
    data::Union{MarkedTimeSeries, Nothing}
    μ::Real
    ϕ::Real
    γ::Real
    δ::Real
    markdens::SupportedMarksDistributions
    ξ::Real
    α::Real
    β::Real
    κ::Real


    function SEMPPExpKern(data::Union{MarkedTimeSeries, Nothing} = nothing; μ::Real = rand(), ϕ::Real = rand(), γ::Real = rand(), δ::Real = rand(), markdens::SupportedMarksDistributions = Distributions.GeneralizedPareto, ξ::Real = rand(), α::Real = rand(), β::Real = rand(), κ::Real = rand())
        if any((μ, ϕ, γ, α, β, δ, κ) .< 0)
            error("paramaters except for ξ must be positive or zero") 
        else 
            new(data, μ, ϕ, γ, δ, markdens, ξ, α, β, κ)
        end
    end

end


params(sepp::SEMPPExpKern) = Dict(:μ => sepp.μ, :ϕ => sepp.ϕ, :γ => sepp.γ, :δ => sepp.δ, :markdens => sepp.markdens, :ξ => sepp.ξ, :α => sepp.α, :β => sepp.β, :κ => sepp.κ)


"""
    lin_coeff_impact(sepp::SEMPPExpKern)::Real

Get the linear coefficient impact of the marks on the rate.
"""
lin_coeff_impact(sepp::SEMPPExpKern)::Real = sepp.δ


"""
    marks_scale_params(sepp::SEMPPExpKern)::Real

Get the paramaters weighing the volatility function into the scale parameter of the mark density: `\\sigma_t = \\beta + \\alpha \\nu(t)`.
"""
marks_scale_params(sepp::SEMPPExpKern)::Real = (sepp.α, sepp.β)


shape(sepp::SEMPPExpKern) = sepp.κ


decay(sepp::SEMPPExpKern) = sepp.ξ