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