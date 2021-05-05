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

    times_plt = plot(x = anytimes, y = points, Guide.yticks(ticks = nothing, label = false), Geom.point)

    rate_plt = plot(x = anytimes, y = lamb, color = [color("black")], Geom.line)

    return vstack(rate_plt, times_plt)
end