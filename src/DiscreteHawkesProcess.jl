"""
volfunc(when, pp, γ, δ)

*volfunc* is the term in the rate of a SEPP that corresponds to the self-excitement. In Li2020 it is denoted by ν.
"""
function volfunc(when::AbstractVector, pp::PP, γ::Real, δ::Real = 0)
    
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


"""
    discrete_negloglik(pp, markdens, μ, ϕ, γ, δ, ξ, β, α)

*discrete_negloglik* calculates the negative log-likelihood of the SE(M)PP with the events and paramaters passed in arguments. 

TODO : if discrete_negloglik is to be exported, it might be useful to add methods so that discrete_negloglik can take a model as an argument instead of its paramaters
"""
function discrete_negloglik(pp::PP;  μ::Real, ϕ::Real, γ::Real)        # one method for point process with or without marks (model without marks)
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

function discrete_negloglik(pp::PP, sepp::SEPPExpKern)
    θ = params(sepp)
    return discrete_negloglik(pp ; θ...)
end


# one methode for marked process

function discrete_negloglik(mpp::MarkedPointProcess, markdens ;  μ::Real, ϕ::Real, γ::Real, δ::Real = 0, ξ::Real, α::Real, β::Real, κ::Real = 1)
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
            return logcdf(markdens(0, sig_mark[1], ξ), sig_mark[2])
        else        # EGPD case
            return logcdf(markdens(sig_mark[1], ξ, κ), sig_mark[2])
        end
    end


    sig_marks = hcat(σ, marks)
    mark_contrib = log_cdf_markdens.(eachrow(sig_marks))

    term2 = sum(mark_contrib)
    
    return -term1 - term2
end


function discrete_negloglik(mpp::MarkedPointProcess, sempp::SEMPPExpKern)
    θ = params(sempp)
    markdens = θ[:markdens]
    delete!(θ, :markdens)
    return discrete_negloglik(mpp, markdens ; θ...)
end


function discrete_fit!(sepp::SEPP, pp::PP) # generic method either to fit a pp whithout marks or to fit the ground process of an mpp
    model = Model(Ipopt.Optimizer)
    θ = params(sepp)

    function to_min(μ, ϕ, γ)
        return discrete_negloglik(pp, μ = μ, ϕ = ϕ, γ = γ)
    end

    JuMP.register(model, :to_min, 3, to_min, autodiff=true)

    @variables(model, begin
        mu >= 0, (start = θ[:μ])
        phi >= 0, (start = θ[:ϕ])
        gamma >= 0, (start = θ[:γ])
    end)
    @NLobjective(model, Min, to_min(mu, phi, gamma))
    
    optimize!(model)

    sepp.μ = value(mu)
    sepp.ϕ = value(phi)
    sepp.γ = value(gamma)

    return objective_value(model)
end


function discrete_fit!(sempp::SEMPPExpKern, mpp::MarkedPointProcess, bounds::Union{Vector{Real}, Nothing} = nothing) # default xi >= 0
    model = Model(Ipopt.Optimizer)
    θ = params(sempp)
    markdens = θ[:markdens]


    function to_min_GPD(μ, ϕ, γ, δ, ξ, α, β)
        return discrete_negloglik(mpp, Distributions.GeneralizedPareto, μ = μ, ϕ = ϕ, γ = γ, δ = δ, ξ = ξ, α = α, β = β)
    end


    function to_min_EGPD(μ, ϕ, γ, δ, ξ, α, β, κ)
        return discrete_negloglik(mpp, markdens, μ = μ, ϕ = ϕ, γ = γ, δ = δ, ξ = ξ, α = α, β = β, κ = κ)
    end


    if markdens == Distributions.GeneralizedPareto
        JuMP.register(model, :to_min_GPD, 7, to_min_GPD, autodiff=true)
    else
        JuMP.register(model, :to_min_EGPD, 8, to_min_EGPD, autodiff=true)
    end


    @variables(model, begin
        mu >= 0, (start = θ[:μ])
        phi >= 0, (start = θ[:ϕ])
        gamma >= 0, (start = θ[:γ])
        delta >= 0, (start = θ[:δ])
        alpha >= 0, (start = θ[:α])
        beta >= 0, (start = θ[:β])
    end)

    if markdens != Distributions.GeneralizedPareto
        @variable(model, kappa >= 0, start = θ[:κ])
    end

    if isnothing(bounds)
        @variable(model, xi >= 0, start = θ[:ξ])
    else
        @variable(model, bounds[1] <= xi <= bounds[2], start = θ[:ξ])
    end


    if markdens == Distributions.GeneralizedPareto
        @NLobjective(model, Min, to_min_GPD(mu, phi, gamma, delta, xi, alpha, beta))
    else
        @NLobjective(model, Min, to_min_EGPD(mu, phi, gamma, delta, xi, alpha, beta, kappa))
    end
    
    optimize!(model)

    sempp.μ = value(mu)
    sempp.ϕ = value(phi)
    sempp.γ = value(gamma)
    sempp.δ = value(gamma)
    sempp.ξ = value(gamma)
    sempp.α = value(gamma)
    sempp.β = value(gamma)

    if markdens != Distributions.GeneralizedPareto
        sempp.κ = value(kappa)
    end

    return objective_value(model)
end