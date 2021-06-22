macro def_discrete_fit_sepp(params_to_fit)
    l = (qn -> qn.value).(params_to_fit.args)
    x = Symbol[]
    variables = :(@variable(model, $(Symbol("y1")) >= 0, start = θ[$(l[1])]))
    
    fitted_params = :(setproperty!(sepp, $(l[1]), value($(Symbol("y1")))))

    h_it = Expr[]

    kw_expr_it = Expr[]

    opt_var = Symbol[]

    for i in 1:length(l)
        push!(x, Symbol("x$i"))
        push!(opt_var, Symbol("y$i"))
        push!(kw_expr_it, :($(l[i]) => $(Symbol("x$i"))))
        push!(h_it, :(getproperty(sepp, $(l[i]))))
    end

    for i in 2:length(l)
        variables = :($variables ; @variable(model, $(Symbol("y$i")) >= 0, start = θ[$(l[i])]))

        fitted_params = :($fitted_params ; $(:(setproperty!(sepp, $(l[i]), value($(Symbol("y$i")))))))
    end

    kw_expr = :(kw = Dict($(kw_expr_it...)))

    to_min_expr = Expr(:function, Expr(:call, :to_min, x...), quote
        default_kw = Dict{Symbol, Real}(:μ => θ[:μ], :ϕ => θ[:ϕ], :γ => θ[:γ])
        $kw_expr
        final_kw = merge(default_kw, kw)
        return discrete_negloglik(ts; final_kw...)
    end)

    h = :([$(h_it...)])

    m = esc.(l)

    return Expr(:function, Expr(:call, :($(esc(:discrete_fit!))), Expr(Symbol("::"), :sepp, :SEPPExpKern), m...),quote
        
        ts = sepp.data
        isnothing(ts) && error("No data in model, can't fit")
        
        model = Model(Ipopt.Optimizer)
        θ = params(sepp)
        
        $to_min_expr
        
        JuMP.register(model, :to_min, length($params_to_fit) , to_min, autodiff=true)

        $variables
        
        @NLobjective(model, Min, to_min($(opt_var...)))
        
        optimize!(model)

        $fitted_params

        to_hess(θ) = to_min(θ...)

        hess = ForwardDiff.hessian(to_hess, $h)

        sepp.cov_mat = inv(hess)

        return objective_value(model)
        
    end)
end


macro def_discrete_fit_sempp(params_to_fit)
    l = (qn -> qn.value).(params_to_fit.args)
    x = Symbol[]
    variables = :(@variable(model, $(Symbol("y1")) >= 0, start = θ[$(l[1])]))
    
    fitted_params = :(setproperty!(sepp, $(l[1]), value($(Symbol("y1")))))

    h_it = Expr[]

    kw_expr_it = Expr[]

    opt_var = Symbol[]

    for i in 1:length(l)
        push!(x, Symbol("x$i"))
        push!(opt_var, Symbol("y$i"))
        push!(kw_expr_it, :($(l[i]) => $(Symbol("x$i"))))
        push!(h_it, :(getproperty(sepp, $(l[i]))))
    end

    for i in 2:length(l)
        variables = :($variables ; @variable(model, $(Symbol("y$i")) >= 0, start = θ[$(l[i])]))

        fitted_params = :($fitted_params ; $(:(setproperty!(sepp, $(l[i]), value($(Symbol("y$i")))))))
    end

    kw_expr = :(kw = Dict($(kw_expr_it...)))

    to_min_expr = Expr(:function, Expr(:call, :to_min, x...), quote
        default_kw = Dict{Symbol, Real}(:μ => θ[:μ], :ϕ => θ[:ϕ], :γ => θ[:γ], :δ => θ[:δ], :ξ => θ[:ξ], :β => θ[:β], :α => θ[:α], :κ => θ[:κ])
        $kw_expr
        final_kw = merge(default_kw, kw)
        return discrete_negloglik(mts, markdens; final_kw...)
    end)

    h = :([$(h_it...)])

    m = esc.(l)

    return Expr(:function, Expr(:call, :($(esc(:discrete_fit!))), Expr(Symbol("::"), :sempp, :SEMPPExpKern), m...),quote
        
        mts = sempp.data
        isnothing(mts) && error("No data in model, can't fit")

        model = Model(Ipopt.Optimizer)
        θ = params(sempp)
        markdens = θ[:markdens]
        
        $to_min_expr
        
        JuMP.register(model, :to_min, length($params_to_fit) , to_min, autodiff=true)

        $variables
        
        @NLobjective(model, Min, to_min($(opt_var...)))
        
        optimize!(model)

        $fitted_params

        to_hess(θ) = to_min(θ...)

        hess = ForwardDiff.hessian(to_hess, $h)

        sempp.cov_mat = inv(hess)

        return objective_value(model)
        
    end)
end