"""
    function negloglik(ts::TS; μ::Real = rand(), ϕ::Real = rand(), γ::Real = rand())::Real

Compute the negative log-likelihood of ts with parameters in kwargs.
"""
function negloglik(ts::TS; μ::Real = rand(), ϕ::Real = rand(), γ::Real = rand())::Real
    tst = [μ, ϕ, γ] .< 0
    any(tst) && (return Inf)

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

#=
function nlL1(mts::MarkedTimeSeries; μ::Real = rand(), ϕ::Real = rand(), γ::Real = rand(), δ::Real = 0)
    tst = [μ, ϕ, γ, δ] .< 0
    w = [:μ, :ϕ, :γ, :δ][tst]
    any(tst) && ((μ, ϕ, γ, δ) = abs.((μ, ϕ, γ, δ)) ; @warn string(string(["$symb " for symb in w]...), "must be positive or zero, taking absolute value"))

    times = first(mts.times) isa TimeType ? Dates.value.(Date.(mts.times)) : mts.times
    starttime = first(times)
    endtime = last(times)
    T = endtime - starttime

    marks = mts.marks

    vol = volfunc(times, mts, γ, δ)
    
    term1 = sum(log.(μ .+ ϕ .* vol))
    term2 = μ * T + ϕ/γ * sum((1 .+ δ .* marks) .* (1 .- exp.(-γ .* (endtime .- times))))
    
    return term2 - term1
end
=#

"""
    negloglik(mts::MarkedTimeSeries, markdens ; μ::Real = rand(), ϕ::Real = rand(), γ::Real = rand(), δ::Real = 0, ξ::Real = rand(), α::Real = rand(), β::Real = rand(), κ::Real = 1)::Real

Compute the negative log-likelihood of the marked time series with parmaters in kwargs.
"""
function negloglik(mts::MarkedTimeSeries, markdens::SupportedMarksDistributions, impact_func::Function = (d -> 1 + tanh(d)) ; μ::Real = rand(), ϕ::Real = rand(), γ::Real = rand(), δ::Real = 0, ξ::Real = rand(), α::Real = rand(), β::Real = rand(), κ::Real = 1)::Real
    tst = [μ, ϕ, γ, δ, β, α, κ] .< 0
    any(tst) && (return Inf)

    times = first(mts.times) isa TimeType ? Dates.value.(DateTime.(mts.times)) ./ (1000*3600*24) : mts.times
    starttime = first(times)
    endtime = last(times)
    T = endtime - starttime

    marks = mts.marks

    vol = volfunc(times, mts, γ, δ, impact_func)
    
    term1 = sum(log.(μ .+ ϕ .* vol))
    term2 = μ * T + ϕ/γ * sum((impact_func(δ .* marks)) .* (1 .- exp.(-γ .* (endtime .- times))))

    σ = β .+ α .* vol
    sig_marks = hcat(σ, marks)

    if markdens == Distributions.GeneralizedPareto
        mark_contrib = (sig_mark -> logpdf(markdens(0, sig_mark[1], ξ), sig_mark[2])).(eachrow(sig_marks))
    else        # EGPD case
        mark_contrib = (sig_mark -> logpdf(markdens(sig_mark[1], ξ, κ), sig_mark[2])).(eachrow(sig_marks))
    end

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
    impact_func = θ[:impact_function]
    delete!(θ, :markdens)
    delete!(θ, :impact_function)
    return negloglik(mts, markdens, impact_func ; θ...)
end


"""
    fit!(sepp::SEPPExpKern)::Real

Fit a self exciting point process model (whithout marks) to the time series in data, based on MLE. 

Note that if the time series has marks, this method ignores them, that is modelling the ground process as independent of the marks.
Returns the minimal Log-likelihood found.
""" 
function fit! end

include("macro_def_fits.jl")

@def_fit_sepp [:s1]
@def_fit_sepp [:s1, :s2]
@def_fit_sepp [:s1, :s2, :s3]


function fit!(sepp::SEPPExpKern)
    return fit!(sepp, :μ, :ϕ, :γ)
end


@def_fit_sempp [:s1]
@def_fit_sempp [:s1, :s2]
@def_fit_sempp [:s1, :s2, :s3]
@def_fit_sempp [:s1, :s2, :s3, :s4]
@def_fit_sempp [:s1, :s2, :s3, :s4, :s5]
@def_fit_sempp [:s1, :s2, :s3, :s4, :s5, :s6]
@def_fit_sempp [:s1, :s2, :s3, :s4, :s5, :s6, :s7]
@def_fit_sempp [:s1, :s2, :s3, :s4, :s5, :s6, :s7, :s8]

function fit!(sempp::SEMPPExpKern)
    return fit!(sempp, :μ, :ϕ, :γ, :δ, :ξ, :α, :β)
end

#=
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
    
    x = [sepp.μ, sepp.ϕ, sepp.γ]

    to_hess(θ) = to_min(θ...)

    hess = ForwardDiff.hessian(to_hess, x)

    sepp.cov_mat = inv(hess)

    return objective_value(model)
end


"""
    fit!(sempp::SEMPPExpKern, bounds::Union{Vector{<:Real}, Nothing} = nothing)::Real

Fit a self exciting marked point process model to the marked time series in data, based on MLE. 

Returns the minimal Log-likelihood found.
""" 
function fit!(sempp::SEMPPExpKern, bounds::Union{Vector{<:Real}, Nothing} = nothing)::Real # default xi >= 0
    
    mts = sempp.data
    isnothing(mts) && error("No data in model, can't fit")

    markdens = sempp.markdens

    if markdens == Distributions.GeneralizedPareto
        df = DataFrame(Marks = mts.marks)
        gml = gpfit(df, :Marks)
        sempp.β = exp(gml.θ̂[1])
        sempp.ξ = gml.θ̂[2]
    else
        egppower = EGPD.EGPpowerfit(mts.marks)
        sempp.ξ = egppower.ξ
        sempp.β = egppower.σ
        sempp.κ = egppower.κ
    end

    model = Model(Ipopt.Optimizer)
    θ = params(sempp)


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

    if markdens == Distributions.GeneralizedPareto
        to_hess1(θ) = to_min_GPD(θ...)
        hess = ForwardDiff.hessian(to_hess1, x)
    else
        sempp.κ = value(kappa)
        push!(x, sempp.κ)
        to_hess2(θ) = to_min_EGPD(θ...)
        hess = ForwardDiff.hessian(to_hess2, x)
    end

    sepp.cov_mat = inv(hess)

    return objective_value(model)
end
=#