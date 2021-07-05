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
        push!(sims, discrete_simulation(sepp; start_time = 0, end_time = horiz))
    end

    warn = 0

    function count(ts::TS)::Real
        times = ts.times
        res = 0
        c = 1

        for i in 2:length(times)
            if c == 0 || times[i-1] == times[i] - 1
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


function monte_carlo_return_period(sepp::SEMPPExpKern; magnitude::Real, r::Integer = 7, horiz::Integer = 100, M::Integer = 500)::Real
    sims = MarkedTimeSeries[]

    for i in 1:M
        push!(sims, discrete_simulation(sepp; start_time = 0, end_time = horiz))
    end

    warn = 0

    function count(mts::MarkedTimeSeries)::Real
        times = mts.times
        marks = mts.marks
        
        res = 0
        c = 0
        marks[1] >= magnitude && (c = 1)

        for i in 2:length(times)
            if (marks[1] >= magnitude) && (c == 0 || times[i-1] == times[i] - 1)
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
        mag_min = min(mag_min, minimum(sim.marks))
        mag_max = max(mag_max, maximum(sim.marks))
        push!(sims, sim)
    end

    println(mag_min)
    println(mag_max)

    function count(mts::MarkedTimeSeries, magnitude::Real)::Real
        times = mts.times
        marks = mts.marks
        
        res = 0
        c = 0
        marks[1] >= magnitude && (c = 1)

        for i in 2:length(times)
            if (marks[1] >= magnitude) && (c == 0 || times[i-1] == times[i] - 1)
                c+=1
            else
                c >= r && (res += 1)
                c=0
            end      
        end

        c >= r && (res += 1)

        res == 0 && (res = 1)

        return horiz/res
    end

    ret_per_min = median((x -> count(x, mag_min)).(sims))
    ret_per_max = median((x -> count(x, mag_max)).(sims))
    
    while ret_per_min != ret_per_max
        mag = 0.5 * (mag_max + mag_min)
        ret_per = median((x -> count(x, mag)).(sims))

        if ret_per > r
            mag_min = mag
            ret_per_min = ret_per
        else
            mag_max = mag    
            ret_per_max = ret_per
        end
    end

    return mag_min, mag_max, ret_per_min
end