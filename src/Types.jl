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


abstract type SEPP end

baseline(sepp::SEPP) = sepp.μ

# les deux noms suivants sont à vérifier

exp_decay(sepp::SEPP) = sepp.γ
self_excitement_factor(sepp::SEPP) = sepp.ϕ

#=
DiscreteSEPPExpKern(μ, ϕ, γ)

*DiscreteSEPPExpKern* is a discrete SEPP with exponential kernel model type with rate λ (for all t, there is an event with probability λ(t)) viz. :

λ(t) = μ + ϕ \sum_{t_k < t} \exp(γ(t-t_k))

where t_k are the timestamps of all events.
=#

mutable struct DiscreteSEPPExpKern <: SEPP
    μ::Real
    ϕ::Real
    γ::Real

    function DiscreteSEPPExpKern(μ::Real = rand(), ϕ::Real = rand(), γ::Real = rand()) 
        if any((μ, ϕ, γ) .< 0) 
            error("paramaters must be positive or zero")
        else 
            new(μ, ϕ, γ)
        end
    end

end

params(sepp::DiscreteSEPPExpKern) = Dict(:μ => sepp.μ, :ϕ => sepp.ϕ, :γ => sepp.γ)

#=
DiscreteSEMPPExpKern(μ, ϕ, γ, markdens, α, β)

*DiscreteSEMPPExpKern* is a discrete SEMPP with exponential kernel model type such that for all t,m , there is an event in t with mark m with probability λ_g(t)f(m|t)) where :

λ_g(t) = μ + ϕ * ν(t)

with ν(t) = \sum_{t_k < t} (1 + δ*m_k) (\exp(γ(t-t_k))

f(m|t) is markdens with scale σ_t = β + α * ν(t)

t_k are the timestamps of all events and m_k their marks.
=#

mutable struct DiscreteSEMPPExpKern <: SEPP
    μ::Real
    ϕ::Real
    γ::Real
    δ::Real
    markdens::Any
    ξ::Real
    α::Real
    β::Real
    κ::Real


    function DiscreteSEMPPExpKern(μ::Real = rand(), ϕ::Real = rand(), γ::Real = rand(), δ::Real = rand(), markdens::Any = Distributions.GeneralizedPareto, ξ::Real = rand(), α::Real = rand(), β::Real = rand(), κ::Real = rand())
        if any((μ, ϕ, γ, α, β, δ, κ) .< 0)
            error("paramaters except for ξ must be positive or zero") 
        else 
            new(μ, ϕ, γ, δ, markdens, ξ, α, β, κ)
        end
    end

end

# noms à vérifier pour les parametres

params(sepp::DiscreteSEMPPExpKern) = Dict(:μ => sepp.μ, :ϕ => sepp.ϕ, :γ => sepp.γ, :δ => sepp.δ, :markdens => sepp.markdens, :ξ => sepp.ξ, :α => sepp.α, :β => sepp.β, :κ => sepp.κ)
lin_coeff_impact(sepp::DiscreteSEMPPExpKern) = sepp.δ
marks_scale_params(sepp::DiscreteSEMPPExpKern) = (sepp.α, sepp.β)
shape(sepp::DiscreteSEMPPExpKern) = sepp.κ
decay(sepp::DiscreteSEMPPExpKern) = sepp.ξ
