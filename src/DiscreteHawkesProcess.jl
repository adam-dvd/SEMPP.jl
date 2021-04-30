

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