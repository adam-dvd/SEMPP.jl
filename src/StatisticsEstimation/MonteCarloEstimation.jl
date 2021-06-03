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


function monte_carlo_return_period(sepp::SEPP, r::Integer = 7, horiz::Integer = 100, M::Integer = 500)::Real
    sims = TS[]

    for i in 1:M
        push!(sims, discrete_simulation(sepp; start_time = 0, end_time = horiz))
    end

    function count(ts::TS)::Real
        times = ts.times
        res = 0
        c = 1

        for i in 2:length(times)
            if c == 0 || times[i-1] == times[i] - 1
                c+=1
            else
                c >= r || (res += 1)
                c=0
            end      
        end

        c >= r || (res += 1)

        return 100/res
    end

    ret_per = median(count.(sims))

    return ret_per
end


function monte_carlo_return_period(sepp::SEMPPExpKern, magnitude::Real, r::Integer = 7, horiz::Integer = 100, M::Integer = 500)::Real
    sims = MarkedTimeSeries[]

    for i in 1:M
        push!(sims, discrete_simulation(sepp; start_time = 0, end_time = horiz))
    end

    function count(mts::MarkedTimeSeries)::Real
        times = mts.times
        marks = mts.marks
        
        res = 0
        c = 0
        marks[1] >= magnitude || (c = 1)

        for i in 2:length(times)
            if (marks[1] >= magnitude) && (c == 0 || times[i-1] == times[i] - 1)
                c+=1
            else
                c >= r || (res += 1)
                c=0
            end      
        end

        c >= r || (res += 1)

        return res/100
    end

    ret_per = median(count.(sims))

    return ret_per
end