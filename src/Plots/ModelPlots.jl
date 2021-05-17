function rate_plot(sepp::SEPP; step = nothing, from_idx = nothing, to_idx = nothing)

    ts = sepp.data
    isnothing(ts) && error("No data in the model, can't plot")

    isnothing(from_idx) && (from_idx = 1)
    isnothing(to_idx) && (to_idx = length(ts))

    times = ts.times
    starttime = times[from_idx]
    endtime = times[to_idx]
    anytimes = isnothing(step) ? (starttime:oneunit(starttime-endtime):endtime) : (starttime:step:endtime)

    μ = sepp.μ
    ϕ = sepp.ϕ
    γ = sepp.γ
    δ = sepp isa SEMPPExpKern ? sepp.δ : 0

    lamb = μ .+ ϕ .* volfunc(anytimes, ts, γ, δ)

    points = fill(0, (to_idx-from_idx,))
    
    rate_layer = layer(x = anytimes, y = lamb, color = [color("black")], Geom.line)

    plt = plot(x = times[from_idx:to_idx], y = points, size = [0mm], Geom.point, Guide.xrug, Theme(grid_line_width=0mm), rate_layer)

    return plt
end


function marks_plot(mts::MarkedTimeSeries)
    marks = mts.marks
    times = mts.times

    plt = plot(x=times, y=marks, Geom.point) 
    
    for k in 1:length(times)
        push!(plt, layer(x=[times[k],times[k]], y=[0, marks[k]], Geom.line))
    end

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