function rate_plot(sepp::SEPP; step = nothing, from = nothing, to = nothing)

    pp = sepp.data
    isnothing(pp) && error("No data in the model, can't plot")

    times = pp.times
    starttime = isnothing(from) ? start_time(pp) : from
    endtime = isnothing(to) ? end_time(pp) : to
    anytimes = isnothing(step) ? (starttime:oneunit(starttime-endtime):endtime) : (starttime:step:endtime)

    μ = sepp.μ
    ϕ = sepp.ϕ
    γ = sepp.γ
    δ = sepp isa SEMPPExpKern ? sepp.δ : 0

    lamb = μ .+ ϕ .* volfunc(anytimes, pp, γ, δ)

    points = fill(0, size(times))

    rate_layer = layer(x = anytimes, y = lamb, color = [color("black")], Geom.line)

    plt = plot(x = times, y = points, size = [0mm], Geom.point, Guide.xrug, Theme(grid_line_width=0mm), rate_layer)

    return plt
end


function marks_plot(mpp::MarkedPointProcess)
    marks = mpp.marks
    times = mpp.times

    plt = plot(x=times, y=marks, Geom.point) 
    
    for k in 1:length(times)
        push!(plt, layer(x=[times[k],times[k]], y=[0, marks[k]], Geom.line))
    end

    return plt
end

function marked_rate_plot(sepp::SEPP, step=nothing)

    mpp = MarkedPointProcess(sepp.data)
    isnothing(mpp) && error("No data in the model, can't plot")
    
    marks = mpp.marks
    times = mpp.times
    starttime = start_time(mpp)
    endtime = end_time(mpp)
    anytimes = isnothing(step) ? (starttime:oneunit(endtime-starttime):endtime) : (starttime:step:endtime)

    μ = sepp.μ
    ϕ = sepp.ϕ
    γ = sepp.γ
    δ = sepp isa SEMPPExpKern ? sepp.δ : 0

    lamb = μ .+ ϕ .* volfunc(anytimes, pp, γ, δ)

    rate_layer = layer(x = anytimes, y = lamb, color = [color("black")], Geom.line)

    plt = plot(x=times, y=marks, Geom.point) 
    
    for k in 1:length(times)
        push!(plt, layer(x=[times[k],times[k]], y=[0, marks[k]], Geom.line))
    end

    push!(plt, rate_layer)

    return plt
end


function marked_rate_plot(sempp::SEMPPExpKern, step= nothing)
    mpp = sempp.data
    isnothing(mpp) && error("No data in the model, can't plot marks")
     
    return marked_rate_plot(sempp, step)
end