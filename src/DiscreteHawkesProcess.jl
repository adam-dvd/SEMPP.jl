
function volfunc(when::AbstractVector, times::AbstractVector, γ::Real, marks::AbstractVector = zero(times), δ::Real = 0)
    
    γ >= 0 && error("γ must be positive or zero")

    (all(mpp .== 0) && δ != 0) && warn("no marks but δ non zero")

    size(times)==size(marks) && error("times and marks must be the same size")

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


function negloglik(times::AbstractVector, marks::AbstractVector, markdens::Function ;  μ::Real, ϕ::Real, γ::Real, δ::Real = 0, ξ::Real, β::Real, α::Real, κ::Real = 1)
    
    (all((μ, ϕ, γ, δ, β, α, κ) .>= 0)) && ((μ, ϕ, γ, δ, β, α, κ) = abs.((μ, ϕ, γ, δ, β, α, κ)) ; warn("all paramaters except for ξ must be positive or zero, taking absolute value"))

    endtime = last(times)
    starttime = first(times)
    anytimes= starttime:endtime

    vol = volfunc(anytimes, times, γ, marks, δ)     # ν function in Li2020
    intens =μ .+ ϕ .* vol       # rate λ in Li2020
    prob = 1 - exp.(-intens )       # probability for an event to happen
    t_idx = findall(in(times), anytimes) 
    prob_1 = prob[t_idx]        # probability of the events that happened to happen
    prob_0 = 1 .- prob[findall(!in(times), anytimes)] # probability of the events that didn't happen not to happen
    term1 = sum(log.(prob_1)) + sum(log.(prob_0))

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