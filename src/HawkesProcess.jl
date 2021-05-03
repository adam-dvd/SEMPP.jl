function negloglik(pp::PointProcess ; μ::Real = rand(), ϕ::Real = rand(), γ::Real = rand())
    
    times = pp.times
    starttime = start_time(pp)
    endtime = end_time(pp)
    T = endtime - starttime

    vol = volfunc(times, pp, γ)
    term1 = sum(log.(μ .+ ϕ .* vol))
    term2 = μ * T + ϕ/γ * sum(1 .- exp.(-γ .* (endtime .- times)))

    return term2 - term1
end


function negloglik(pp::PointProcess, sepp::SEPPExpKern)
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
            return logcdf(markdens(0, ξ, sig_mark[1]), sig_mark[2])
        else        # EGPD case
            return logcdf(markdens(sig_mark[1], ξ, κ), sig_mark[2])
        end
    end


    sig_marks = hcat(σ, marks)
    mark_contrib = log_cdf_markdens.(eachrow(sig_marks))

    term3 = sum(mark_contrib)

    return term2 - term1 - term3

end

function discrete_negloglik(mpp::MarkedPointProcess, sempp::SEMPPExpKern)
    θ = params(sempp)
    markdens = θ[:markdens]
    delete!(θ, :markdens)
    return negloglik(mpp, markdens ; θ...)
end


function fit!(sepp::SEPP, pp::PointProcess)
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

