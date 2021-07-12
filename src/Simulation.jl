"""
    simulation::TS

Simulate event records from the model.
"""
function simulation end


function simulation(sepp::SEPPExpKern; start_time::Real = 0, end_time::Real=1000, history_time_series::Union{TimeSeries, Nothing} = nothing)::TimeSeries
    
    first(history_time_series.times) isa TimeType && (error("TimeType are discrete, try discrete_simulation or use a continuous time series."))

    μ = sepp.μ
    ϕ = sepp.ϕ
    γ = sepp.γ
    
    last_h = -Inf   # keeps track of the last time stamp that is not simulated

    if isnothing(history_time_series)
        times = Float64[]
        ts = TimeSeries(times)
        real_start = start_time
        real_end = end_time
    else
        times = copy(history_time_series.times)
        last_h = last(times)
        real_start = last_h + start_time
        real_end = last_h + end_time 
        ts = TimeSeries(times)
    end

    t = real_start
    λ = μ + ϕ * volfunc([t], ts, γ)[1]

    while t <= real_end
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

    ts = TimeSeries(times[times .> last_h])

    return ts
end


function simulation(sempp::SEMPPExpKern; start_time::Real = 0, end_time::Real=1000, history_marked_time_series::Union{MarkedTimeSeries, Nothing} = nothing)::MarkedTimeSeries
    
    first(history_marked_time_series.times) isa TimeType && (error("TimeType are discrete, try discrete_simulation or use a continuous time series."))

    μ = sempp.μ
    ϕ = sempp.ϕ
    γ = sempp.γ
    δ = sempp.δ
    markdens = sempp.markdens
    ξ = sempp.ξ
    α = sempp.α
    β = sempp.β
    κ = sempp.κ

    last_h = -Inf   # keeps track of the last time stamp that is not simulated

    if isnothing(history_marked_time_series)
        times = Float64[]
        marks = Float64[]
        mts = MarkedTimeSeries(times, marks)
    else
        times = copy(history_marked_time_series.times)
        marks = copy(history_marked_time_series.marks)
        last_h = last(times)
        real_start = last_h + start_time
        real_end = last_h + end_time 
        mts = MarkedTimeSeries(times, marks)
    end

    t = real_start
    λ = μ + ϕ * volfunc([t], mts, γ, δ)[1]

    while t <= real_end
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

    idx = times .> last_h
    mts = MarkedTimeSeries(times[idx], marks[idx])

    return mts
end


"""
    discrete_simulation::TS

Simulate event records from the model on a discrete timeline.
"""
function discrete_simulation end


function discrete_simulation(sepp::SEPPExpKern; start_time::Integer = 1, end_time::Integer = 1000, history_time_series::Union{TimeSeries, Nothing} = nothing)::TimeSeries
    μ = sepp.μ
    ϕ = sepp.ϕ
    γ = sepp.γ

    time_bool = false
    last_h = -Inf   # keeps track of the last time stamp that is not simulated

    if isnothing(history_time_series)
        times = Integer[]
        ts = TimeSeries(times)
        real_start = start_time
        real_end = end_time
    else
        time_bool = first(history_time_series.times) isa TimeType
        times = time_bool ? Dates.value.(Date.(history_time_series.times)) : copy(history_time_series.times)
        last_h = last(times)
        real_start = last_h + start_time
        real_end = last_h + end_time
        ts = TimeSeries(times)
    end
    
    for t in real_start:real_end
        λ = μ + ϕ * volfunc([t], ts, γ)[1]
        prob = 1 - exp(-λ)

        if rand() <= prob
            push!(times, t)
            ts = TimeSeries(times)
        end
    end

    ts = TimeSeries(times[times .> last_h])

    time_bool && (ts = TimeSeries(Date.(times)))

    return ts
end


function discrete_simulation(sempp::SEMPPExpKern; start_time::Int = 1, end_time::Int=1000, history_marked_time_series::Union{MarkedTimeSeries, Nothing} = nothing)::MarkedTimeSeries
    μ = sempp.μ
    ϕ = sempp.ϕ
    γ = sempp.γ
    δ = sempp.δ
    markdens = sempp.markdens
    ξ = sempp.ξ
    α = sempp.α
    β = sempp.β
    κ = sempp.κ

    time_bool = false
    last_h = -Inf   # keeps track of the last time stamp that is not simulated

    if isnothing(history_marked_time_series)
        times = Float64[]
        marks = Float64[]
        mts = MarkedTimeSeries(times, marks)
    else
        time_bool = first(history_marked_time_series.times) isa TimeType
        times = time_bool ? Dates.value.(Date.(history_marked_time_series.times)) : copy(history_marked_time_series.times)
        marks = copy(history_marked_time_series.marks)
        last_h = last(times)
        real_start = last_h + start_time
        real_end = last_h + end_time
        mts = MarkedTimeSeries(times, marks)
    end
    
    for t in real_start:real_end
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

    idx = times .> last_h
    mts = MarkedTimeSeries(times[idx], marks[idx])

    time_bool && (mts = MarkedTimeSeries(Date.(times), marks))

    return mts
end