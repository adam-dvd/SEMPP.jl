"""
    simulation::TS

Simulate event records from the model.
"""
function simulation end


function simulation(sepp::SEPPExpKern; start_time::Real = 0, end_time::Real=1000, history_time_series::Union{TimeSeries, Nothing} = nothing)::TimeSeries
    
    μ = sepp.μ
    ϕ = sepp.ϕ
    γ = sepp.γ

    first(history_time_series.times) isa TimeType && (error("TimeType are discrete, try discrete_simulation or use a continuous time series."))

    if isnothing(history_time_series)
        times = Float64[]
        ts = TimeSeries(times)
    else
        times = history_time_series.times
        start_time = last(times) 
        end_time = start_time + end_time 
        ts = TimeSeries(times)
    end

    t = start_time
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


function simulation(sempp::SEMPPExpKern; start_time::Real = 0, end_time::Real=1000, history_marked_time_series::Union{MarkedTimeSeries, Nothing} = nothing)::MarkedTimeSeries
    
    μ = sempp.μ
    ϕ = sempp.ϕ
    γ = sempp.γ
    δ = sempp.δ
    markdens = sempp.markdens
    ξ = sempp.ξ
    α = sempp.α
    β = sempp.β
    κ = sempp.κ

    first(history_marked_time_series.times) isa TimeType && (error("TimeType are discrete, try discrete_simulation or use a continuous time series."))

    if isnothing(history_marked_time_series)
        times = Float64[]
        marks = Float64[]
        mts = MarkedTimeSeries(times, marks)
    else
        times = history_marked_time_series.times
        marks = history_marked_time_series.marks
        start_time = last(times) 
        end_time = start_time + end_time 
        mts = MarkedTimeSeries(times, marks)
    end

    t = start_time
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


function discrete_simulation(sepp::SEPPExpKern; start_time::Int = 0, end_time::Int=1000, history_time_series::Union{TimeSeries, Nothing} = nothing)::TimeSeries
    μ = sepp.μ
    ϕ = sepp.ϕ
    γ = sepp.γ

    time_bool = first(history_time_series.times) isa TimeType

    if isnothing(history_time_series)
        times = Float64[]
        ts = TimeSeries(times)
    else
        times = time_bool ? Dates.value.(Date.(history_time_series.times)) : history_time_series.times
        start_time = last(times) + 1
        end_time = start_time + end_time - 1
        ts = TimeSeries(times)
    end
    
    for t in start_time:end_time
        λ = μ + ϕ * volfunc([t], ts, γ)[1]
        prob = 1 - exp(-λ)

        if rand() <= prob
            push!(times, t)
            ts = TimeSeries(times)
        end
    end

    time_bool && ts = TimeSeries(Date.(times))

    return ts
end


function discrete_simulation(sempp::SEMPPExpKern; start_time::Int = 0, end_time::Int=1000, history_marked_time_series::Union{MarkedTimeSeries, Nothing} = nothing)::MarkedTimeSeries
    μ = sempp.μ
    ϕ = sempp.ϕ
    γ = sempp.γ
    δ = sempp.δ
    markdens = sempp.markdens
    ξ = sempp.ξ
    α = sempp.α
    β = sempp.β
    κ = sempp.κ

    time_bool = first(history_marked_time_series.times) isa TimeType

    if isnothing(history_marked_time_series)
        times = Float64[]
        marks = Float64[]
        mts = MarkedTimeSeries(times, marks)
    else
        times = time_bool ? Dates.value.(Date.(history_marked_time_series.times)) : history_marked_time_series.times
        marks = history_marked_time_series.marks
        start_time = last(times) + 1
        end_time = start_time + end_time - 1
        mts = MarkedTimeSeries(times, marks)
    end
    
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

    time_bool && mts = MarkedTimeSeries(Date.(times), marks)

    return mts
end