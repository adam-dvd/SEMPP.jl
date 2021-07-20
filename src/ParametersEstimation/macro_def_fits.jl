macro def_fit_sepp(params_to_fit)
    l = (qn -> qn.value).(params_to_fit.args)
    n = length(l)
    x0_expr = Expr[]
    fitted_params = :(setproperty!(sepp, $(l[1]), Optim.minimizer(res)[1]))

    h_it = Expr[]

    kw_expr_it = Expr[]

    for i in 1:n
        push!(x0_expr, :(θ[$(l[i])]))
        push!(kw_expr_it, :($(l[i]) => x[$i]))
        push!(h_it, :(getproperty(sepp, $(l[i]))))
    end

    for i in 2:n
        fitted_params = :($fitted_params ; $(:(setproperty!(sepp, $(l[i]), Optim.minimizer(res)[$i]))))
    end

    kw_expr = :(kw = Dict($(kw_expr_it...)))

    to_min_expr = Expr(:function, Expr(:call, :to_min, :x), quote
        default_kw = Dict{Symbol, Real}(:μ => θ[:μ], :ϕ => θ[:ϕ], :γ => θ[:γ])
        $kw_expr
        final_kw = merge(default_kw, kw)
        return negloglik(ts; final_kw...)
    end)

    h = :([$(h_it...)])

    x0_expr = :([$(x0_expr...)])

    m = esc.(l)

    return Expr(:function, Expr(:call, :($(esc(:fit!))), Expr(Symbol("::"), :sepp, :SEPPExpKern), m...),quote
        
        ts = sepp.data
        isnothing(ts) && error("No data in model, can't fit")

        θ = params(sepp)
        
        $to_min_expr

        lower = zeros($n)
        upper = fill(Inf, $n)
        
        res = optimize(to_min, lower, upper, $x0_expr)

        $fitted_params

        hess = ForwardDiff.hessian(to_min, $h)

        sepp.cov_mat = inv(hess)

        return res
        
    end)
end


macro def_fit_sempp(params_to_fit)
    l = (qn -> qn.value).(params_to_fit.args)
    n = length(l)
    x0_expr = Expr[]
    fitted_params = :(setproperty!(sempp, $(l[1]), Optim.minimizer(res)[1]))

    h_it = Expr[]

    kw_expr_it = Expr[]

    for i in 1:n
        push!(x0_expr, :(θ[$(l[i])]))
        push!(kw_expr_it, :($(l[i]) => x[$i]))
        push!(h_it, :(getproperty(sempp, $(l[i]))))
    end

    for i in 2:n
        fitted_params = :($fitted_params ; $(:(setproperty!(sempp, $(l[i]), Optim.minimizer(res)[$i]))))
    end

    kw_expr = :(kw = Dict($(kw_expr_it...)))

    to_min_expr = Expr(:function, Expr(:call, :to_min, :x), quote
        default_kw = Dict{Symbol, Real}(:μ => θ[:μ], :ϕ => θ[:ϕ], :γ => θ[:γ])
        $kw_expr
        final_kw = merge(default_kw, kw)
        return negloglik(mts, markdens, impact_func; final_kw...)
    end)

    h = :([$(h_it...)])

    x0_expr = :([$(x0_expr...)])

    m = esc.(l)

    return Expr(:function, Expr(:call, :($(esc(:fit!))), Expr(:parameters, Expr(:kw, :lower, zeros(n)), Expr(:kw, :upper, fill(Inf, n))), Expr(Symbol("::"), :sempp, :SEMPPExpKern), m...),quote
        
        mts = sempp.data
        isnothing(mts) && error("No data in model, can't fit")

        θ = params(sempp)
        markdens = θ[:markdens]
        impact_func = θ[:impact_function]
        
        $to_min_expr
        
        res = optimize(to_min, lower, upper, $x0_expr)

        $fitted_params

        hess = ForwardDiff.hessian(to_min, $h)

        sempp.cov_mat = inv(hess)

        return res
        
    end)
end