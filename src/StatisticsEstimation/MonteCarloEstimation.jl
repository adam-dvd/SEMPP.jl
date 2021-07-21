#=
function q(len::Int = 7, ynb::Int = 2, horiz::Int = 100, M::Int = 500)
    
end
=#


"""
    monte_carlo_return_period(sepp::SEPP, r::Integer = 7, horiz::Integer = 100, M::Integer = 500)::Real

Compute the retrun period of an event of *r* days by a monte calro method.

The computation is based on *M* simulation to the horizon *horiz*. 
If the model is an SEMPP and magnitude is given by the user than the method only consider events above *magnitude*.
"""
function monte_carlo_return_period end


function monte_carlo_return_period(sepp::SEPP; r::Integer = 7, horiz::Integer = 100, M::Integer = 500)::Real
    sims = TS[]

    for i in 1:M
        push!(sims, discrete_simulation(sepp; start_time = 1, end_time = horiz))
    end

    warn = 0

    function count(ts::TS)::Real
        times = copy(ts.times)
        res = 0
        c = 1

        for i in 2:length(times)
            if c == 0 || times[i-1] == times[i] - 1
                c += 1
            else
                c >= r && (res += 1)
                c = 0
            end      
        end

        c >= r && (res += 1)

        res == 0 && (res = 1 ; warn += 1)

        return horiz/res
    end

    ret_per = median(count.(sims))

    warn == 0 || @warn "no events within horizon for $warn simulation(s)"

    return ret_per
end


function monte_carlo_return_period(sempp::SEMPPExpKern; magnitude::Real, r::Integer = 7, horiz::Integer = 100, M::Integer = 500)::Real
    sims = MarkedTimeSeries[]

    for i in 1:M
        push!(sims, discrete_simulation(sempp; start_time = 1, end_time = horiz))
    end

    warn = 0

    function count(mts::MarkedTimeSeries)::Real
        times = copy(mts.times)
        marks = copy(mts.marks)
        
        res = 0
        c = 0
        marks[1] >= magnitude && (c = 1)

        for i in 2:length(times)
            if (marks[i] >= magnitude) && (c == 0 || times[i-1] == times[i] - 1)
                c+=1
            else
                c >= r && (res += 1)
                c=0
            end      
        end

        c >= r && (res += 1)

        res == 0 && (res = 1 ; warn += 1)

        return horiz/res
    end

    ret_per = median(count.(sims))

    warn == 0 || @warn "no events within horizon for $warn simulation(s)"

    return ret_per
end


function rQy(sempp::SEMPPExpKern, r::Integer = 7, y::Real = 2; horiz::Union{Real, Nothing} = nothing, M::Integer = 500)
    sims = MarkedTimeSeries[]
    horiz = isnothing(horiz) ? 5*y*365 : Int(horiz*365)

    mag_min = Inf
    mag_max = -Inf

    for i in 1:M
        sim = discrete_simulation(sempp; start_time = 0, end_time = horiz)
        mag_min = minimum(vcat(sim.marks, mag_min))
        mag_max = maximum(vcat(sim.marks, mag_max))
        push!(sims, sim)
    end

    function count(mts::MarkedTimeSeries, magnitude::Real)::Real
        times = copy(mts.times)
        marks = copy(mts.marks)
        
        res = 0
        c = 0

        if !isempty(times)
            marks[1] >= magnitude && (c = 1)

            for i in 2:length(times)
                if (marks[i] >= magnitude) && (c == 0 || times[i-1] == times[i] - 1)
                    c+=1
                else
                    c >= r && (res += 1)
                    c=0
                end      
            end
        end
        
        c >= r && (res += 1)

        res == 0 && (res = 1)

        return horiz/res
    end

    ret_per_min = median((x -> count(x, mag_min)).(sims))
    ret_per_max = median((x -> count(x, mag_max)).(sims))

    println(ret_per_min)
    println(ret_per_max)
    
    while abs(ret_per_max - ret_per_min) >= 1
        mag = 0.5 * (mag_max + mag_min)
        ret_per = median((x -> count(x, mag)).(sims))

        if ret_per > y*365
            mag_max = mag    
            ret_per_max = ret_per
        else
            mag_min = mag
            ret_per_min = ret_per
        end
    end

    return 0.5*(mag_min + mag_max), ret_per_min
end


function discrete_forecast(sepp::SEPPExpKern; start_time::Integer = 1, end_time::Integer = 1000, history::Union{TimeSeries, Nothing} = nothing, M::Integer = 500)
    sims = TimeSeries[]

    if !isnothing(history)
        (first(history.times) isa TimeType) && (history = TimeSeries(Dates.value.(Date.(history.times))))
    end

    last_h = last(history.times)

    for i in 1:M
        push!(sims, discrete_simulation(sepp, start_time = start_time, end_time = end_time, history_time_series = history))
    end

    p = zeros(Float64, end_time - start_time + 1)

    real_start = last_h + start_time

    for ts in sims
        for t in ts.times
            p[t - real_start + 1] += 1
        end
    end

    p = p ./ M

    return p
end


function discrete_forecast(sempp::SEMPPExpKern; start_time::Integer = 1, end_time::Integer = 1000, history::Union{MarkedTimeSeries, Nothing} = nothing, magnitudes::AbstractVector{<:Real} = [], M::Integer = 500)
    sims = MarkedTimeSeries[]

    if !isnothing(history)
        (first(history.times) isa TimeType) && (history = MarkedTimeSeries(Dates.value.(Date.(history.times)), history.marks))
    end

    for i in 1:M
        push!(sims, discrete_simulation(sempp, start_time = start_time, end_time = end_time, history_time_series = history))
    end

    p = zeros(Float64, end_time - start_time + 1, length(magnitudes) - 1)

    for mts in sims
        extended_marks = fill(-Inf, end_time - start_time + 1)
        extended_marks[mts.times .- start_time .+ 1] = mts.marks
        
        for i in 1:(end_time - start_time + 1)
            m_idx = sum(magnitudes .<= extended_marks[i])

            if (m_idx > 0) 
                p[i, m_idx] += 1
            else
                p[i, 1] += 1
            end
        end
    end

    p = p ./ M

    return p
end