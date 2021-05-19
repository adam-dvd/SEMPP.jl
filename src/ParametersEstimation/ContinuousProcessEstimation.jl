"""
    function negloglik(ts::TS; μ::Real = rand(), ϕ::Real = rand(), γ::Real = rand())::Real

Compute the negative log-likelihood of ts with parameters in kwargs.
"""
function negloglik(ts::TS; μ::Real = rand(), ϕ::Real = rand(), γ::Real = rand())::Real
    tst = [μ, ϕ, γ] .< 0
    w = [:μ, :ϕ, :γ][tst]
    any(tst) && ((μ, ϕ, γ) = abs.((μ, ϕ, γ)) ; @warn string(string(["$symb " for symb in w]...), "must be positive or zero, taking absolute value"))
    
    times = first(ts.times) isa TimeType ? Dates.value.(DateTime.(ts.times)) ./ (1000*3600*24) : ts.times
    starttime = first(times)
    endtime = last(times)
    T = endtime - starttime

    vol = volfunc(times, ts, γ)
    term1 = sum(log.(μ .+ ϕ .* vol))
    term2 = μ * T + ϕ/γ * sum(1 .- exp.(-γ .* (endtime .- times)))

    return term2 - term1
end


"""
    negloglik(sepp::SEPPExpKern)::Real

Compute the negative log-likelihood of the time series in data with parmaters of the sepp model.
"""
function negloglik(sepp::SEPPExpKern)::Real
    ts = sepp.data
    isnothing(ts) && error("No data in model, can't compute log-likelihood")

    θ = params(sepp)
    return negloglik(ts ; θ...)
end


"""
    negloglik(mts::MarkedTimeSeries, markdens ; μ::Real = rand(), ϕ::Real = rand(), γ::Real = rand(), δ::Real = 0, ξ::Real = rand(), α::Real = rand(), β::Real = rand(), κ::Real = 1)::Real

Compute the negative log-likelihood of the marked time series with parmaters in kwargs.
"""
function negloglik(mts::MarkedTimeSeries, markdens ; μ::Real = rand(), ϕ::Real = rand(), γ::Real = rand(), δ::Real = 0, ξ::Real = rand(), α::Real = rand(), β::Real = rand(), κ::Real = 1)::Real
    tst = [μ, ϕ, γ, δ, β, α, κ] .< 0
    w = [:μ, :ϕ, :γ, :δ, :β, :α, :κ][tst]
    any(tst) && ((μ, ϕ, γ, δ, β, α, κ) = abs.((μ, ϕ, γ, δ, β, α, κ)) ; @warn string(string(["$symb " for symb in w]...), "must be positive or zero, taking absolute value"))

    times = first(mts.times) isa TimeType ? Dates.value.(DateTime.(mts.times)) ./ (1000*3600*24) : mts.times
    starttime = first(times)
    endtime = last(times)
    T = endtime - starttime

    marks = mts.marks

    vol = volfunc(times, mts, γ)
    
    term1 = sum(log.(μ .+ ϕ .* vol))
    term2 = μ * T + ϕ/γ * sum((1 .+ δ .* marks) .* (1 .- exp.(-γ .* (endtime .- times))))

    σ = β .+ α .* vol


    function log_cdf_markdens(sig_mark)
        if markdens == Distributions.GeneralizedPareto
            return logcdf(markdens(0, sig_mark[1], ξ), sig_mark[2])
        else        # EGPD case
            return logcdf(markdens(sig_mark[1], ξ, κ), sig_mark[2])
        end
    end


    sig_marks = hcat(σ, marks)
    mark_contrib = log_cdf_markdens.(eachrow(sig_marks))

    term3 = sum(mark_contrib)

    return term2 - term1 - term3

end


"""
    negloglik(sempp::SEMPPExpKern)::Real

 Compute the negative log-likelihood of the marked time series in data with parmaters of the sempp model.
"""
function negloglik(sempp::SEMPPExpKern)::Real
    mts = sempp.data
    isnothing(mts) && error("No data in model, can't compute log-likelihood")

    θ = params(sempp)
    markdens = θ[:markdens]
    delete!(θ, :markdens)
    return negloglik(mts, markdens ; θ...)
end


"""
    fit!(sepp::SEPPExpKern)::Real

Fit a self exciting point process model (whithout marks) to the time series in data, based on MLE. 

Note that if the time series has marks, this method ignores them, that is modelling the ground process as independent of the marks.
Returns the minimal Log-likelihood found.
""" 
function fit!(sepp::SEPPExpKern)::Real
    ts = sepp.data
    isnothing(ts) && error("No data in model, can't fit")

    model = Model(Ipopt.Optimizer)
    θ = params(sepp)

    function to_min(μ, ϕ, γ)
        return negloglik(ts, μ = μ, ϕ = ϕ, γ = γ)
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

#=
"""
    fit!(sempp::SEMPPExpKern, bounds::Union{Vector{<:Real}, Nothing} = nothing)::Real

Fit a self exciting marked point process model to the marked time series in data, based on MLE. 

Returns the minimal Log-likelihood found.
""" 
function fit!(sempp::SEMPPExpKern, bounds::Union{Vector{<:Real}, Nothing} = nothing)::Real # default xi >= 0
    mts = sempp.data
    isnothing(mts) && error("No data in model, can't fit")

    model = Model(Ipopt.Optimizer)
    θ = params(sempp)
    markdens = θ[:markdens]


    function to_min_GPD(μ, ϕ, γ, δ, ξ, α, β)
        return negloglik(mts, Distributions.GeneralizedPareto, μ = μ, ϕ = ϕ, γ = γ, δ = δ, ξ = ξ, α = α, β = β)
    end


    function to_min_EGPD(μ, ϕ, γ, δ, ξ, α, β, κ)
        return negloglik(mts, markdens, μ = μ, ϕ = ϕ, γ = γ, δ = δ, ξ = ξ, α = α, β = β, κ = κ)
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
=#

function fit!(sempp::SEMPPExpKern, bounds::Union{Vector{<:Real}, Nothing} = nothing)::Real # default xi >= 0
    mts = sempp.data
    isnothing(mts) && error("No data in model, can't fit")

    sepp = SEPPExpKern(TimeSeries(copy(mts.times)))
    fit!(sepp)
    sempp.μ = sepp.μ
    sempp.ϕ = sepp.ϕ
    sempp.γ = sepp.γ

    if markdens == Distributions.GeneralizedPareto
        df = DataFrame(Marks = mts.marks)
        gml = gpfit(df, :Marks)
        sempp.ξ, sempp.α = gml.θ̂
    else
        egppower = EGPD.EGPpowerfit(mts.marks)
        sempp.ξ = egppower.ξ
        sempp.α = egppower.σ
        sempp.κ = egppower.κ
    end

    model = Model(Ipopt.Optimizer)
    θ = params(sempp)
    markdens = θ[:markdens]


    function to_min_GPD(μ, ϕ, γ, δ, ξ, α, β)
        return negloglik(mts, Distributions.GeneralizedPareto, μ = μ, ϕ = ϕ, γ = γ, δ = δ, ξ = ξ, α = α, β = β)
    end


    function to_min_EGPD(μ, ϕ, γ, δ, ξ, α, β, κ)
        return negloglik(mts, markdens, μ = μ, ϕ = ϕ, γ = γ, δ = δ, ξ = ξ, α = α, β = β, κ = κ)
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