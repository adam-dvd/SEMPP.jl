function negloglik(pp::PP; μ::Real = rand(), ϕ::Real = rand(), γ::Real = rand())
    
    times = pp.times
    starttime = start_time(pp)
    endtime = end_time(pp)
    T = endtime - starttime

    vol = volfunc(times, pp, γ)
    term1 = sum(log.(μ .+ ϕ .* vol))
    term2 = μ * T + ϕ/γ * sum(1 .- exp.(-γ .* (endtime .- times)))

    return term2 - term1
end


function negloglik(sepp::SEPPExpKern)
    pp = sepp.data
    isnothing(pp) && error("No data in model, can't compute log-likelihood")

    θ = params(sepp)
    return negloglik(pp ; θ...)
end


function negloglik(mpp::MarkedPointProcess, markdens ; μ::Real = rand(), ϕ::Real = rand(), γ::Real = rand(), δ::Real = 0, ξ::Real = rand(), α::Real = rand(), β::Real = rand(), κ::Real = 1)
    
    times = mpp.times
    starttime = start_time(mpp)
    endtime = end_time(mpp)
    T = endtime - starttime

    marks = mpp.marks

    vol = volfunc(times, mpp, γ)
    
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

function negloglik(sempp::SEMPPExpKern)
    mpp = sempp.data
    isnothing(mpp) && error("No data in model, can't compute log-likelihood")

    θ = params(sempp)
    markdens = θ[:markdens]
    delete!(θ, :markdens)
    return negloglik(mpp, markdens ; θ...)
end


function fit!(sepp::SEPPExpKern)
    pp = sepp.data
    isnothing(pp) && error("No data in model, can't fit")

    model = Model(Ipopt.Optimizer)
    θ = params(sepp)

    function to_min(μ, ϕ, γ)
        return negloglik(pp, μ = μ, ϕ = ϕ, γ = γ)
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


function fit!(sempp::SEMPPExpKern, bounds::Union{Vector{<:Real}, Nothing} = nothing) # default xi >= 0
    mpp = sempp.data
    isnothing(mpp) && error("No data in model, can't fit")

    model = Model(Ipopt.Optimizer)
    θ = params(sempp)
    markdens = θ[:markdens]


    function to_min_GPD(μ, ϕ, γ, δ, ξ, α, β)
        return negloglik(mpp, Distributions.GeneralizedPareto, μ, ϕ, γ, δ, ξ, α, β)
    end


    function to_min_EGPD(μ, ϕ, γ, δ, ξ, α, β, κ)
        return negloglik(mpp, markdens, μ, ϕ, γ, δ, ξ, α, β, κ)
    end


    if markdens == Distributions.GeneralizedPareto
        JuMP.register(model, :to_min_GPD, 7, to_min, autodiff=true)
    else
        JuMP.register(model, :to_min_EGPD, 8, to_min, autodiff=true)
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

