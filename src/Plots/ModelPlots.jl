function rate_plot(sepp::SEPP; step = nothing, from_idx = nothing, to_idx = nothing)

    ts = sepp.data
    isnothing(ts) && error("No data in the model, can't plot")

    times = ts.times
    isnothing(from_idx) && (from_idx = 1)
    isnothing(to_idx) && (to_idx = length(times))
    starttime = times[from_idx]
    endtime = times[to_idx]

    ((starttime isa TimeType) && !(step isa Union{TimeType, Nothing})) && (@warn "If times are TimeType step must be TimeType, ignoring step value."; step=nothing)

    if isnothing(step)
        if starttime isa TimeType
            anytimes = DateTime(starttime):Dates.Hour(1):DateTime(endtime)
        else
            anytimes = (starttime:oneunit(starttime-endtime):endtime)
        end
    else
        anytimes = (starttime:step:endtime)
    end

    μ = sepp.μ
    ϕ = sepp.ϕ
    γ = sepp.γ
    δ = sepp isa SEMPPExpKern ? sepp.δ : 0

    lamb = μ .+ ϕ .* volfunc(anytimes, ts, γ, δ)

    points = fill(0, (to_idx-from_idx + 1,))
    
    rate_layer = layer(x = anytimes, y = lamb, color = [color("black")], Geom.line)

    points_layer = layer(x = times[from_idx:to_idx], y = points, alpha = [0.7], Geom.point)

    plt = plot(rate_layer, points_layer, Guide.xlabel("time"), Guide.ylabel("intensity"), Theme(grid_line_width=0mm, highlight_width = 0mm))

    return plt
end


function marks_plot(mts::MarkedTimeSeries; from_idx::Union{Nothing, Integer} = nothing, to_idx::Union{Nothing, Integer} = nothing)
    marks = mts.marks
    times = mts.times

    isnothing(from_idx) && (from_idx = 1)
    isnothing(to_idx) && (to_idx = length(times))

    plt = plot(x=times[from_idx:to_idx], y=marks[from_idx:to_idx], Geom.hair, Geom.point, Guide.xlabel("time"), Guide.ylabel("marks")) 

    return plt
end


function marks_plot(sempp::SEMPPExpKern)
    mts = sempp.data
    isnothing(mts) && error("No data in the model, can't plot")
    
    return marks_plot(mts)
end


function marked_rate_plot(sempp::SEMPPExpKern; step = nothing, from_idx = nothing, to_idx = nothing)

    mts = sempp.data
    isnothing(mts) && error("No data in the model, can't plot")

    rate_plt = rate_plot(sempp, step = step, from_idx = from_idx, to_idx = to_idx)

    marks_plt = marks_plot(mts, from_idx = from_idx, to_idx = to_idx)

    plt = vstack(rate_plt, marks_plt) 

    return plt
end