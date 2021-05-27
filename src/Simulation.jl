"""
    simulation::TS

Simulate event records from the model.
"""
function simulation end


function simulation(sepp::SEPPExpKern, end_time::Real=1000)::TimeSeries
    
    μ = sepp.μ
    ϕ = sepp.ϕ
    γ = sepp.γ

    times = Float64[]
    ts = TimeSeries(times)
    t = 0
    λ = μ + ϕ * volfunc([t], ts, γ)[1]

    while t <= end_time
        M = λ
        s = rand(Distributions.Exponential(1/λ))
        t = t + s
        λ = μ + ϕ * volfunc([t], ts, γ)[1]
        p = rand()

        if p <= λ/M 
            push!(times, t)
            ts = TimeSeries(times)
        end

    end

    return ts
end


function simulation(sempp::SEMPPExpKern, end_time::Real=10)::MarkedTimeSeries
    
    μ = sempp.μ
    ϕ = sempp.ϕ
    γ = sempp.γ
    δ = sempp.δ
    markdens = sempp.markdens
    ξ = sempp.ξ
    α = sempp.α
    β = sempp.β
    κ = sempp.κ

    times = Float64[]
    marks = Float64[]
    mts = MarkedTimeSeries(times, marks)
    t = 0
    λ = μ + ϕ * volfunc([t], mts, γ, δ)[1]

    while t <= end_time
        M = λ
        s = rand(Distributions.Exponential(1/λ))
        t = t + s
        vol = volfunc([t], mts, γ, δ)
        λ = μ + ϕ * vol[1]
        p = rand()

        if p <= λ/M 
            push!(times, t)
            σ = β + α * vol[1]

            if markdens == Distributions.GeneralizedPareto
                m = rand(Distributions.GeneralizedPareto(0, σ, ξ))
            else
                m = rand(EGPD.EGPpower(σ, ξ, κ))
            end

            push!(marks, m)
            mts = MarkedTimeSeries(times, marks)
        end

    end

    return mts
end


"""
    simulation::TS

Simulate event records from the model on a discrete timeline.
"""
function discrete_simulation end


function discrete_simulation(sepp::SEPPExpKern; start_time::Int = 0, end_time::Int=100)::TimeSeries
    μ = sepp.μ
    ϕ = sepp.ϕ
    γ = sepp.γ

    times = Float64[]
    ts = TimeSeries(times)
    
    for t in start_time:end_time
        λ = μ + ϕ * volfunc([t], ts, γ)[1]
        prob = 1 - exp(-λ)

        if rand() <= prob
            push!(times, t)
            ts = TimeSeries(times)
        end
    end

    return ts
end


function discrete_simulation(sempp::SEMPPExpKern; start_time::Int = 0, end_time::Int=100)::MarkedTimeSeries
    μ = sempp.μ
    ϕ = sempp.ϕ
    γ = sempp.γ
    δ = sempp.δ
    markdens = sempp.markdens
    ξ = sempp.ξ
    α = sempp.α
    β = sempp.β
    κ = sempp.κ

    times = Float64[]
    marks = Float64[]
    mts = MarkedTimeSeries(times, marks)
    
    for t in start_time:end_time
        vol = volfunc([t], mts, γ, δ)[1]
        λ = μ + ϕ * vol
        prob = 1 - exp(-λ)

        if rand() <= prob
            push!(times, t)

            σ = β + α * vol[1]

            if markdens == Distributions.GeneralizedPareto
                m = rand(Distributions.GeneralizedPareto(0, σ, ξ))
            else
                m = rand(EGPD.EGPpower(σ, ξ, κ))
            end

            push!(marks, m)
            
            mts = MarkedTimeSeries(times, marks)
        end
    end

    return mts
end