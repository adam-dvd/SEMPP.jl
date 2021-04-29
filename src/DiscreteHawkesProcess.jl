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


#=
volfunc(when, pp, γ, δ)

*volfunc* is the term in the rate of a SEPP that corresponds to the self-excitement. In Li2020 it is denoted by ν.
=#

function volfunc(when::AbstractVector, pp::SEMPP.PP, γ::Real, δ::Real = 0)
    
    γ < 0 && error("γ must be positive or zero")

    (pp isa SEMPP.PointProcess && δ != 0) && (@warn "no marks but δ non zero")

    times = pp.times
    marks = pp isa SEMPP.MarkedPointProcess ? pp.marks : zero(times)

    mpp = hcat(times, marks)


    function self_ex(t)

        mpp_to_t = mpp[times .< t, :]

        function term(markedpoint)
            return (1+ δ*markedpoint[2])*exp(-γ*(t-markedpoint[1]))
        end


        return sum(term.(eachrow(mpp_to_t)))
    end


    return self_ex.(when)
end

#=
negloglik(pp, markdens, μ, ϕ, γ, δ, ξ, β, α)

*negloglik* calculates the negative log-likelihood of the SE(M)PP with the events and paramaters passed in arguments. 

TODO : if negloglik is to be exported, it might be useful to add methods so that negloglik can take a model as an argument instead of its paramaters
=#

# one method for point process without marks

function negloglik(pp::PointProcess;  μ::Real, ϕ::Real, γ::Real)
    tst = [μ, ϕ, γ] .< 0
    w = [:μ, :ϕ, :γ][tst]
    any(tst) && ((μ, ϕ, γ) = abs.((μ, ϕ, γ)) ; @warn string(string(["$symb " for symb in w]...), "must be positive or zero, taking absolute value"))

    times = pp.times
    endtime = end_time(pp)
    starttime = start_time(pp)
    anytimes= starttime:endtime

    vol = volfunc(anytimes, pp, γ)     # ν function in Li2020
    intens = μ .+ ϕ .* vol       # rate λ in Li2020
    prob = 1 .- exp.(-intens )       # probability for an event to happen
    t_idx = findall(in(times), anytimes) 
    prob_1 = prob[t_idx]        # probability of the events that happened to happen
    prob_0 = 1 .- prob[findall(!in(times), anytimes)] # probability of the events that didn't happen not to happen
    term1 = sum(log.(prob_1)) + sum(log.(prob_0))
 
    return (-term1)
end 

function negloglik(pp::PointProcess, sepp::DiscreteSEPPExpKern)
    θ = params(sepp)
    return negloglik(pp ; θ...)
end


# one methode for marked process

function negloglik(mpp::MarkedPointProcess, markdens ;  μ::Real, ϕ::Real, γ::Real, δ::Real = 0, ξ::Real, α::Real, β::Real, κ::Real = 1)
    tst = [μ, ϕ, γ, δ, β, α, κ] .< 0
    w = [:μ, :ϕ, :γ, :δ, :β, :α, :κ][tst]
    any(tst) && ((μ, ϕ, γ, δ, β, α, κ) = abs.((μ, ϕ, γ, δ, β, α, κ)) ; @warn string(string(["$symb " for symb in w]...), "must be positive or zero, taking absolute value"))

    times = mpp.times
    endtime = end_time(mpp)
    starttime = start_time(mpp)
    anytimes= starttime:endtime

    vol = volfunc(anytimes, mpp, γ, δ)     # ν function in Li2020
    intens =μ .+ ϕ .* vol       # rate λ in Li2020
    prob = 1 .- exp.(-intens )       # probability for an event to happen
    t_idx = findall(in(times), anytimes) 
    prob_1 = prob[t_idx]        # probability of the events that happened to happen
    prob_0 = 1 .- prob[findall(!in(times), anytimes)] # probability of the events that didn't happen not to happen
    term1 = sum(log.(prob_1)) + sum(log.(prob_0))

    σ = β .+ α .* vol[t_idx]
    marks = mpp.marks


    function log_cdf_markdens(sig_mark)
        if markdens == Distributions.GeneralizedPareto
            return logcdf(markdens(0, ξ, sig_mark[1]), sig_mark[2])
        else        # EGPD case
            return logcdf(markdens(sig_mark[1], ξ, κ), sig_mark[2])
        end
    end


    sig_marks = hcat(σ, marks)
    mark_contrib = log_cdf_markdens.(eachrow(sig_marks))

    term2 = sum(mark_contrib)
    
    return -term1 - term2
end


function negloglik(mpp::MarkedPointProcess, sempp::DiscreteSEMPPExpKern)
    θ = params(sempp)
    markdens = θ[:markdens]
    delete!(θ, :markdens)
    return negloglik(mpp, markdens ; θ...)
end


function fit!(sepp::SEPP, pp::PP)
    model = Model(GLPK.Optimizer)
    θ = params(sepp)
    @variables(model, begin
        mu >= 0, (start = θ[:μ])
        phi >= 0, (start = θ[:ϕ])
        gamma >= 0, (start = θ[:γ])
    end)
    @NLobjective(model, Min, negloglik(pp, μ = mu, ϕ = phi, γ = gamma))
    
    optimize!(model)

    sepp.μ = value(μ)
    sepp.ϕ = value(ϕ)
    sepp.γ = value(γ)
    
    println(objective_value(model))

    return
end