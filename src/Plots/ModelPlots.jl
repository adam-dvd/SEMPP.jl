function rate_plot(sepp::SEPP, pp::PP, step::Real = 1)

    times = pp.times
    starttime = start_time(pp)
    endtime = end_time(pp)
    anytimes = starttime:step:endtime

    μ = sepp.μ
    ϕ = sepp.ϕ
    γ = sepp.γ
    δ = sepp isa SEMPPExpKern ? sepp.δ : 0

    lamb = μ .+ ϕ .* volfunc(anytimes, pp, γ, δ)

    points = zero(times)

    rate_layer = layer(x = anytimes, y = lamb, color = [color("black")], Geom.line)

    plt = plot(x = anytimes, y = points, size = [0mm], Geom.point, Guide.xrug, Theme(grid_line_width=0mm), rate_layer)

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

function marked_rate_plot(sepp::SEPP, mpp::MarkedPointProcess, step::Real = 1)
    
    marks = mpp.marks
    times = mpp.times
    starttime = start_time(mpp)
    endtime = end_time(mpp)
    anytimes = starttime:step:endtime

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