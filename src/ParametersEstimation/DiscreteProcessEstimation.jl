"""
    discrete_negloglik::Real

Compute the negative log-likelihood of the process in argument. 
"""
function discrete_negloglik(ts::TS;  μ::Real, ϕ::Real, γ::Real)::Real        # one method for point process with or without marks (model without marks)
    tst = [μ, ϕ, γ] .< 0
    w = [:μ, :ϕ, :γ][tst]
    any(tst) && ((μ, ϕ, γ) = abs.((μ, ϕ, γ)) ; @warn string(string(["$symb " for symb in w]...), "must be positive or zero, taking absolute value"))

    times = ts.times
    endtime = end_time(ts)
    starttime = start_time(ts)
    anytimes= starttime:oneunit(starttime-endtime):endtime

    vol = volfunc(anytimes, ts, γ)     # ν function in Li2020
    intens = μ .+ ϕ .* vol       # rate λ in Li2020
    prob = 1 .- exp.(-intens )       # probability for an event to happen
    t_idx = findall(in(times), anytimes) 
    prob_1 = prob[t_idx]        # probability of the events that happened to happen
    prob_0 = 1 .- prob[findall(!in(times), anytimes)] # probability of the events that didn't happen not to happen
    term1 = sum(log.(prob_1)) + sum(log.(prob_0))
 
    return (-term1)
end


function discrete_negloglik(sepp::SEPPExpKern)
    ts = sepp.data
    isnothing(ts) && error("No data in model, can't compute log-likelihood")

    θ = params(sepp)
    return discrete_negloglik(ts ; θ...)
end


function discrete_negloglik(mts::MarkedTimeSeries, markdens::SupportedMarksDistributions ;  μ::Real, ϕ::Real, γ::Real, δ::Real = 0, ξ::Real, α::Real, β::Real, κ::Real = 1)
    tst = [μ, ϕ, γ, δ, β, α, κ] .< 0
    w = [:μ, :ϕ, :γ, :δ, :β, :α, :κ][tst]
    any(tst) && ((μ, ϕ, γ, δ, β, α, κ) = abs.((μ, ϕ, γ, δ, β, α, κ)) ; @warn string(string(["$symb " for symb in w]...), "must be positive or zero, taking absolute value"))

    times = mts.times
    endtime = end_time(mts)
    starttime = start_time(mts)
    anytimes= starttime:oneunit(starttime-endtime):endtime

    vol = volfunc(anytimes, mts, γ, δ)     # ν function in Li2020
    intens =μ .+ ϕ .* vol       # rate λ in Li2020
    prob = 1 .- exp.(-intens )       # probability for an event to happen
    t_idx = findall(in(times), anytimes) 
    prob_1 = prob[t_idx]        # probability of the events that happened to happen
    prob_0 = 1 .- prob[findall(!in(times), anytimes)] # probability of the events that didn't happen not to happen
    term1 = sum(log.(prob_1)) + sum(log.(prob_0))

    σ = β .+ α .* vol[t_idx]
    marks = mts.marks

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
    mts = sempp.data
    isnothing(mts) && error("No data in model, can't compute log-likelihood")

    θ = params(sempp)
    markdens = θ[:markdens]
    delete!(θ, :markdens)
    return discrete_negloglik(mts, markdens ; θ...)
end


"""
    discrete_fit!(sepp::SEPP)

Fit a discrete self exciting point process model (whithout marks) to the time series in data, based on MLE. 

Note that if the process has marks, this method ignores them, that is modelling the ground process as independent of the marks.
Returns the minimal Log-likelihood found. 
"""
function discrete_fit!(sepp::SEPP) # generic method either to fit a ts whithout marks or to fit the ground process of an mts
    ts=sepp.data
    isnothing(ts) && error("No data in model, can't fit")

    model = Model(Ipopt.Optimizer)
    θ = params(sepp)

    function to_min(μ, ϕ, γ)
        return discrete_negloglik(ts, μ = μ, ϕ = ϕ, γ = γ)
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
    x = [sepp.μ, sepp.ϕ, sepp.γ]

    d = NLPEvaluator(model)
    MathOptInterface.initialize(d, [:Hess])
    hess_structure = MathOptInterface.hessian_lagrangian_structure(d, values)
    hess_values = zero(hess_structure)
    MathOptInterface.eval_hessian_lagrangian(d, hess_values, x, 1, zero(x))

    sepp.cov_mat = zeros(3, 3)
    
    for i in 1:length(hess_structure)
        sepp.cov_mat[hess_structure[i]] += hess_values[i]
    end

    sepp.cov_mat = inv(sepp.cov_mat)

    return objective_value(model)
end


"""
    discrete_fit!(sempp::SEMPPExpKern, bounds::Union{Vector{<:Real}, Nothing} = nothing)

Fit a discrete self exciting marked point process model with marks distribution either GPD or EGPpower. Based on MLE.
"""
function discrete_fit!(sempp::SEMPPExpKern, bounds::Union{Vector{<:Real}, Nothing} = nothing) # default xi >= 0
    mts = sempp.data
    isnothing(mts) && error("No data in model, can't fit")

    model = Model(Ipopt.Optimizer)
    θ = params(sempp)
    markdens = θ[:markdens]


    function to_min_GPD(μ, ϕ, γ, δ, ξ, α, β)
        return discrete_negloglik(mts, Distributions.GeneralizedPareto, μ = μ, ϕ = ϕ, γ = γ, δ = δ, ξ = ξ, α = α, β = β)
    end


    function to_min_EGPD(μ, ϕ, γ, δ, ξ, α, β, κ)
        return discrete_negloglik(mts, EGPD.EGPpower, μ = μ, ϕ = ϕ, γ = γ, δ = δ, ξ = ξ, α = α, β = β, κ = κ)
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
    sempp.δ = value(delta)
    sempp.ξ = value(xi)
    sempp.α = value(alpha)
    sempp.β = value(beta)
    x = [sempp.μ, sempp.ϕ, sempp.γ, sempp.δ, sempp.ξ, sempp.α, sempp.β]

    if markdens == EGPD.EGPpower
        sempp.κ = value(kappa)
        push!(x, sempp.κ)
    end

    d = NLPEvaluator(model)
    MathOptInterface.initialize(d, [:Hess])
    hess_structure = MathOptInterface.hessian_lagrangian_structure(d, values)
    hess_values = zero(hess_structure)
    MathOptInterface.eval_hessian_lagrangian(d, hess_values, x, 1, zero(x))

    if markdens == Distributions.GeneralizedPareto
        sempp.cov_mat = zeros(7, 7)
    else
        sempp.cov_mat = zeros(8, 8)
    end
    
    for i in 1:length(hess_structure)
        sempp.cov_mat[hess_structure[i]] += hess_values[i]
    end

    sempp.cov_mat = inv(sempp.cov_mat)

    return objective_value(model)
end