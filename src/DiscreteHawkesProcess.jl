
"""
DiscreteSEPPExpKern(μ, ϕ, γ)

*DiscreteSEPPExpKern* is a discrete SEPP with exponential kernel model type with rate λ (for all t, there is an event with probability λ(t)) viz. :

λ(t) = μ + ϕ \sum_{t_k < t} \exp(γ(t-t_k))

where t_k are the timestamps of all events.
"""

struct DiscreteSEPPExpKern
    μ::Real
    ϕ::Real
    γ::Real
    DiscreteSEPPExpKern(μ, ϕ, γ) = any((μ, ϕ, γ) .< 0) ? error("paramaters must be positive or zero") : new(μ, ϕ, γ)
end

"""
DiscreteSEMPPExpKern(μ, ϕ, γ, markdens, α, β)

*DiscreteSEMPPExpKern* is a discrete SEMPP with exponential kernel model type such that for all t,m , there is an event in t with mark m with probability λ_g(t)f(m|t)) where :

λ_g(t) = μ + ϕ * ν(t)

with ν(t) = \sum_{t_k < t} (1 + δ*m_k) (\exp(γ(t-t_k))

f(m|t) is markdens with scale σ_t = β + α * ν(t)

t_k are the timestamps of all events and m_k their marks.
"""

struct DiscreteSEMPPExpKern
    μ::Real
    ϕ::Real
    γ::Real
    markdens::ContinuousUnivariateDistribution
    α::Real
    β::Real
    DiscreteSEPPExpKern(μ, ϕ, γ, markdens, α, β) = any((μ, ϕ, γ, α, β) .< 0) ? error("paramaters must be positive or zero") : new(μ, ϕ, γ, markdens, α, β)
end


"""
volfunc(when, pp, γ, δ)

*volfunc* is the term in the rate of a SEPP that corresponds to the self-excitement. In Li2020 it is denoted by ν.
"""

function volfunc(when::AbstractVector, pp::PP, γ::Real, δ::Real = 0)
    
    γ >= 0 && error("γ must be positive or zero")

    (pp isa PointProcess && δ != 0) && warn("no marks but δ non zero")

    times = pp.times
    marks = pp isa MarkedPointProcess ? pp.marks : zeros(times)

    mpp = hcat(times, marks)

    function self_ex(t)

        mpp_to_t = mpp[times .< t]

        function term(markedpoint)
            return (1+ δ*markedpoint[2])*exp(-γ*(t-markedpoint[1]))
        end


        return sum(term.(eachline(mpp_to_t)))
    end

    return self_ex.(when)
end

"""
negloglik(pp, markdens, μ, ϕ, γ, δ, ξ, β, α)

*negloglik* calculates the negative log-likelihood of the SE(M)PP with the events and paramaters passed in arguments. 

TODO : if negloglik is to be exported, it might be useful to add methods so that negloglik can take a model as an argument instead of its paramaters
"""


function negloglik(pp::PointProcess, markdens::ContinuousUnivariateDistribution ;  μ::Real, ϕ::Real, γ::Real, δ::Real = 0, ξ::Real, β::Real, α::Real, κ::Real = 1)
    
    (all((μ, ϕ, γ, δ, β, α, κ) .>= 0)) && ((μ, ϕ, γ, δ, β, α, κ) = abs.((μ, ϕ, γ, δ, β, α, κ)) ; warn("all paramaters except for ξ must be positive or zero, taking absolute value"))

    endtime = end_time(pp)
    starttime = start_time(pp)
    anytimes= starttime:endtime

    vol = volfunc(anytimes, pp, γ, δ)     # ν function in Li2020
    intens =μ .+ ϕ .* vol       # rate λ in Li2020
    prob = 1 - exp.(-intens )       # probability for an event to happen
    t_idx = findall(in(times), anytimes) 
    prob_1 = prob[t_idx]        # probability of the events that happened to happen
    prob_0 = 1 .- prob[findall(!in(times), anytimes)] # probability of the events that didn't happen not to happen
    term1 = sum(log.(prob_1)) + sum(log.(prob_0))

    if pp isa PointProcess 
        return (-term1)
    end 

    σ = β .+ α .* vol[t_idx]

    if markdens == Distributions.GeneralizedPareto      # GPD case

        function log_cdf_markdens(sig_mark)
            return logcdf(markdens(0, ξ, sig_mark[1]), sig_mark[2])
        end

    else        # EGPD case

        function log_cdf_markdens(sig_mark)
            return logcdf(markdens(sig_mark[1], ξ, κ), sig_mark[2])
        end

    end

    sig_marks = hcat(σ, marks)
    mark_contrib = log_cdf_markdens.(sig_marks)

    term2 = sum(mark_contrib)
    
    return -term1 - term2
end 