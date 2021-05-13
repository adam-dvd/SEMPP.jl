"""
A type for point process might be useful and precise things, but it might conversely impede clarity... not sure to keep it
"""
abstract type PP end


struct PointProcess <: PP
    times::AbstractVector

    function PointProcess(times)
        times = reshape(times, (length(times),))
        return new(times)
    end
end


struct MarkedPointProcess <: PP
    times::AbstractVector
    marks::AbstractVector

    function MarkedPointProcess(times, marks) 
        length(times) == length(marks) || error("times and marks must be the same size")
        times = reshape(times, (length(times),))
        marks = reshape(marks, (length(marks),))
        return new(times, marks)
    end
end

# a method to see a point process as marked with zero marks
function MarkedPointProcess(pp::PP)
    if pp isa PointProcess
        times = pp.times
        marks = fill(0, size(times))
        return MarkedPointProcess(times, marks)
    end
    
    return pp
end


function ground_process(pp::PP)
    pp = PointProcess(pp.times)
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

"""
SEPPExpKern(μ, ϕ, γ)

*SEPPExpKern* is a SEPP with exponential kernel model type with rate λ (for all t, there is an event with probability λ(t)) viz. :
''λ(t) = μ + ϕ \\sum_{t_k < t} \\exp(γ(t-t_k))''
where t_k are the timestamps of all events.
"""
mutable struct SEPPExpKern <: SEPP
    data::Union{PP, Nothing}
    μ::Real
    ϕ::Real
    γ::Real

    function SEPPExpKern(data::Union{PP, Nothing} = nothing; μ::Real = rand(), ϕ::Real = rand(), γ::Real = rand()) 
        if any((μ, ϕ, γ) .< 0) 
            error("paramaters must be positive or zero")
        else 
            new(data,μ, ϕ, γ)
        end
    end

end

params(sepp::SEPPExpKern) = Dict(:μ => sepp.μ, :ϕ => sepp.ϕ, :γ => sepp.γ)


SupportedMarksDistributions = Union{Type{Distributions.GeneralizedPareto}, Type{EGPD.EGPpower}}


"""
SEMPPExpKern(μ, ϕ, γ, markdens, α, β)

*SEMPPExpKern* is a SEMPP with exponential kernel model type such that for all t,m , there is an event in t with mark m with probability λ_g(t)f(m|t)) where :

λ_g(t) = μ + ϕ * ν(t)

with ''ν(t) = \\sum_{t_k < t} (1 + δ*m_k) (\\exp(γ(t-t_k))''

f(m|t) is markdens with scale σ_t = β + α * ν(t)

t_k are the timestamps of all events and m_k their marks.
"""
mutable struct SEMPPExpKern <: SEPP
    data::Union{MarkedPointProcess, Nothing}
    μ::Real
    ϕ::Real
    γ::Real
    δ::Real
    markdens::SupportedMarksDistributions
    ξ::Real
    α::Real
    β::Real
    κ::Real


    function SEMPPExpKern(data::Union{MarkedPointProcess, Nothing} = nothing; μ::Real = rand(), ϕ::Real = rand(), γ::Real = rand(), δ::Real = rand(), markdens::SupportedMarksDistributions = Distributions.GeneralizedPareto, ξ::Real = rand(), α::Real = rand(), β::Real = rand(), κ::Real = rand())
        if any((μ, ϕ, γ, α, β, δ, κ) .< 0)
            error("paramaters except for ξ must be positive or zero") 
        else 
            new(data, μ, ϕ, γ, δ, markdens, ξ, α, β, κ)
        end
    end

end

# noms à vérifier pour les parametres

params(sepp::SEMPPExpKern) = Dict(:μ => sepp.μ, :ϕ => sepp.ϕ, :γ => sepp.γ, :δ => sepp.δ, :markdens => sepp.markdens, :ξ => sepp.ξ, :α => sepp.α, :β => sepp.β, :κ => sepp.κ)
lin_coeff_impact(sepp::SEMPPExpKern) = sepp.δ
marks_scale_params(sepp::SEMPPExpKern) = (sepp.α, sepp.β)
shape(sepp::SEMPPExpKern) = sepp.κ
decay(sepp::SEMPPExpKern) = sepp.ξ


