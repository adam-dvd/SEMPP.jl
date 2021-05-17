function rate_plot(sepp::SEPP; step = nothing, from_idx = nothing, to_idx = nothing)

    ts = sepp.data
    isnothing(ts) && error("No data in the model, can't plot")

    times = ts.times
    isnothing(from_idx) && (from_idx = 1)
    isnothing(to_idx) && (to_idx = length(times))
    starttime = times[from_idx]
    endtime = times[to_idx]
    anytimes = isnothing(step) ? (starttime:oneunit(starttime-endtime):endtime) : (starttime:step:endtime)

    μ = sepp.μ
    ϕ = sepp.ϕ
    γ = sepp.γ
    δ = sepp isa SEMPPExpKern ? sepp.δ : 0

    lamb = μ .+ ϕ .* volfunc(anytimes, ts, γ, δ)

    points = fill(-0.2, (to_idx-from_idx + 1,))
    
    rate_layer = layer(x = anytimes, y = lamb, color = [color("black")], Geom.line)

    points_layer = layer(x = times[from_idx:to_idx], y = points, alpha = [0.8], Geom.point)

    plt = plot(rate_layer, points_layer, Theme(grid_line_width=0mm, highlight_width = 0mm))

    return plt
end


function marks_plot(mts::MarkedTimeSeries)
    marks = mts.marks
    times = mts.times

    plt = plot(x=times, y=marks, Geom.hair, Geom.point) 

    return plt
end


function marked_rate_plot(sepp::SEPP, step=nothing)

    mts = MarkedTimeSeries(sepp.data)
    isnothing(mts) && error("No data in the model, can't plot")
    
    marks = mts.marks
    times = mts.times
    starttime = start_time(mts)
    endtime = end_time(mts)
    anytimes = isnothing(step) ? (starttime:oneunit(endtime-starttime):endtime) : (starttime:step:endtime)

    μ = sepp.μ
    ϕ = sepp.ϕ
    γ = sepp.γ
    δ = sepp isa SEMPPExpKern ? sepp.δ : 0

    lamb = μ .+ ϕ .* volfunc(anytimes, ts, γ, δ)

    rate_layer = layer(x = anytimes, y = lamb, color = [color("black")], Geom.line)

    plt = plot(x=times, y=marks, Geom.point) 
    
    for k in 1:length(times)
        push!(plt, layer(x=[times[k],times[k]], y=[0, marks[k]], Geom.line))
    end

    push!(plt, rate_layer)

    return plt
end


function marked_rate_plot(sempp::SEMPPExpKern, step= nothing)
    mts = sempp.data
    isnothing(mts) && error("No data in the model, can't plot marks")
     
    return marked_rate_plot(sempp, step)
end