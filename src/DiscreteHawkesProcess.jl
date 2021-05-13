"""
    volfunc(when, pp, γ, δ)

compute the term in the rate of a SEPP that corresponds to the self-excitement, denoted by ν in Li2020.
"""
function volfunc(when::AbstractVector, pp::PP, γ::Real, δ::Real = 0)
    
    γ < 0 && error("γ must be positive or zero")

    (pp isa SEMPP.PointProcess && δ != 0) && (@warn "no marks but δ non zero")

    first(when) isa TimeType && when = Dates.value.(when)
    times = first(pp.times) isa TimeType ? Dates.value.(pp.times) : pp.times
    marks = pp isa SEMPP.MarkedPointProcess ? pp.marks : fill(0, size(times))

    mpp = hcat(times, marks)


    function self_ex(t)

        mpp_to_t = mpp[times .< t, :]

        term(markedpoint) = (1+ δ*markedpoint[2])*exp(-γ*(t-markedpoint[1]))

        return sum(term.(eachrow(mpp_to_t)))
    end


    return self_ex.(when)
end


"""
    discrete_negloglik(pp, markdens, μ, ϕ, γ, δ, ξ, β, α)

compute the negative log-likelihood of the process in argument. 

TODO : if discrete_negloglik is to be exported, it might be useful to add methods so that discrete_negloglik can take a model as an argument instead of its paramaters
"""
function discrete_negloglik(pp::PP;  μ::Real, ϕ::Real, γ::Real)        # one method for point process with or without marks (model without marks)
    tst = [μ, ϕ, γ] .< 0
    w = [:μ, :ϕ, :γ][tst]
    any(tst) && ((μ, ϕ, γ) = abs.((μ, ϕ, γ)) ; @warn string(string(["$symb " for symb in w]...), "must be positive or zero, taking absolute value"))

    times = pp.times
    endtime = end_time(pp)
    starttime = start_time(pp)
    anytimes= starttime:oneunit(starttime-endtime):endtime

    vol = volfunc(anytimes, pp, γ)     # ν function in Li2020
    intens = μ .+ ϕ .* vol       # rate λ in Li2020
    prob = 1 .- exp.(-intens )       # probability for an event to happen
    t_idx = findall(in(times), anytimes) 
    prob_1 = prob[t_idx]        # probability of the events that happened to happen
    prob_0 = 1 .- prob[findall(!in(times), anytimes)] # probability of the events that didn't happen not to happen
    term1 = sum(log.(prob_1)) + sum(log.(prob_0))
 
    return (-term1)
end 

function discrete_negloglik(sepp::SEPPExpKern)
    pp = sepp.data
    isnothing(pp) && error("No data in model, can't compute log-likelihood")

    θ = params(sepp)
    return discrete_negloglik(pp ; θ...)
end


# one methode for marked process

function discrete_negloglik(mpp::MarkedPointProcess, markdens::SupportedMarksDistributions ;  μ::Real, ϕ::Real, γ::Real, δ::Real = 0, ξ::Real, α::Real, β::Real, κ::Real = 1)
    tst = [μ, ϕ, γ, δ, β, α, κ] .< 0
    w = [:μ, :ϕ, :γ, :δ, :β, :α, :κ][tst]
    any(tst) && ((μ, ϕ, γ, δ, β, α, κ) = abs.((μ, ϕ, γ, δ, β, α, κ)) ; @warn string(string(["$symb " for symb in w]...), "must be positive or zero, taking absolute value"))

    times = mpp.times
    endtime = end_time(mpp)
    starttime = start_time(mpp)
    anytimes= starttime:oneunit(starttime-endtime):endtime

    vol = volfunc(anytimes, mpp, γ, δ)     # ν function in Li2020
    intens =μ .+ ϕ .* vol       # rate λ in Li2020
    prob = 1 .- exp.(-intens )       # probability for an event to happen
    t_idx = findall(in(times), anytimes) 
    prob_1 = prob[t_idx]        # probability of the events that happened to happen
    prob_0 = 1 .- prob[findall(!in(times), anytimes)] # probability of the events that didn't happen not to happen
    term1 = sum(log.(prob_1)) + sum(log.(prob_0))

    σ = β .+ α .* vol[t_idx]
    marks = mpp.marks

    sig_marks = hcat(σ, marks)

    if markdens == Distributions.GeneralizedPareto
        mark_contrib = (sig_mark -> logcdf(markdens(0, sig_mark[1], ξ), sig_mark[2])).(eachrow(sig_marks))
    else        # EGPpower case
        mark_contrib = (sig_mark -> logcdf(markdens(sig_mark[1], ξ, κ), sig_mark[2])).(eachrow(sig_marks))
    end

    term2 = sum(mark_contrib)
    
    return -term1 - term2
end


function discrete_negloglik(sempp::SEMPPExpKern)
    mpp = sempp.data
    isnothing(mpp) && error("No data in model, can't compute log-likelihood")

    θ = params(sempp)
    markdens = θ[:markdens]
    delete!(θ, :markdens)
    return discrete_negloglik(mpp, markdens ; θ...)
end

"""
    discrete_fit!(sepp, pp)

fit a self exciting point process model (whithout marks) to a time series. 

Note that if the process has marks, this method ignores them, that is modelling the ground process as independent of the marks. 
"""
function discrete_fit!(sepp::SEPP) # generic method either to fit a pp whithout marks or to fit the ground process of an mpp
    pp=sepp.data
    isnothing(pp) && error("No data in model, can't fit")

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

"""
    discrete_fit!(sepp, pp)

fit a self exciting marked point process model with marks distribution either GPD or EGPpower.
"""
function discrete_fit!(sempp::SEMPPExpKern, bounds::Union{Vector{<:Real}, Nothing} = nothing) # default xi >= 0
    mpp = sempp.data
    isnothing(mpp) && error("No data in model, can't fit")

    model = Model(Ipopt.Optimizer)
    θ = params(sempp)
    markdens = θ[:markdens]


    function to_min_GPD(μ, ϕ, γ, δ, ξ, α, β)
        return discrete_negloglik(mpp, Distributions.GeneralizedPareto, μ = μ, ϕ = ϕ, γ = γ, δ = δ, ξ = ξ, α = α, β = β)
    end


    function to_min_EGPD(μ, ϕ, γ, δ, ξ, α, β, κ)
        return discrete_negloglik(mpp, EGPD.EGPpower, μ = μ, ϕ = ϕ, γ = γ, δ = δ, ξ = ξ, α = α, β = β, κ = κ)
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

    if markdens == EGPD.EGPpower
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