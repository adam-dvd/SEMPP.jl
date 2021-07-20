macro def_discrete_fit_sepp(params_to_fit)
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
        return discrete_negloglik(ts; final_kw...)
    end)

    h = :([$(h_it...)])

    x0_expr = :([$(x0_expr...)])

    m = esc.(l)

    return Expr(:function, Expr(:call, :($(esc(:discrete_fit!))), Expr(Symbol("::"), :sepp, :SEPPExpKern), m...),quote
        
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

#=
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

=#

macro def_discrete_fit_sempp(params_to_fit)
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
        return discrete_negloglik(mts, markdens, impact_func; final_kw...)
    end)

    h = :([$(h_it...)])

    x0_expr = :([$(x0_expr...)])

    m = esc.(l)

    return Expr(:function, Expr(:call, :($(esc(:discrete_fit!))), Expr(:parameters, Expr(:kw, :lower, zeros(n)), Expr(:kw, :upper, fill(Inf, n))), Expr(Symbol("::"), :sempp, :SEMPPExpKern), m...),quote
        
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

#=
macro def_discrete_fit_sempp(params_to_fit)
    l = (qn -> qn.value).(params_to_fit.args)
    x = Symbol[]
    variables = :(@variable(model, $(Symbol("y1")) >= 0, start = θ[$(l[1])]))
    constraints = :(@constraint(model, $(Symbol("y1")) >= 0))
    
    fitted_params = :(setproperty!(sempp, $(l[1]), value($(Symbol("y1")))))

    h_it = Expr[]

    kw_expr_it = Expr[]

    opt_var = Symbol[]

    for i in 1:length(l)
        push!(x, Symbol("x$i"))
        push!(opt_var, Symbol("y$i"))
        push!(kw_expr_it, :($(l[i]) => $(Symbol("x$i"))))
        push!(h_it, :(getproperty(sempp, $(l[i]))))
    end

    for i in 2:length(l)
        variables = :($variables ; @variable(model, $(Symbol("y$i")) >= 0, start = θ[$(l[i])]))

        constraints = :($constraints ;  @constraint(model, $(Symbol("y$i")) >= 0))

        fitted_params = :($fitted_params ; $(:(setproperty!(sempp, $(l[i]), value($(Symbol("y$i")))))))
    end

    kw_expr = :(kw = Dict($(kw_expr_it...)))

    to_min_expr = Expr(:function, Expr(:call, :to_min, x...), quote
        default_kw = Dict{Symbol, Real}(:μ => θ[:μ], :ϕ => θ[:ϕ], :γ => θ[:γ], :δ => θ[:δ], :ξ => θ[:ξ], :β => θ[:β], :α => θ[:α], :κ => θ[:κ])
        $kw_expr
        final_kw = merge(default_kw, kw)
        return discrete_negloglik(mts, markdens, impact_func; final_kw...)
    end)

    h = :([$(h_it...)])

    m = esc.(l)

    return Expr(:function, Expr(:call, :($(esc(:discrete_fit!))), Expr(:parameters, Expr(:kw, :tol, 10^(-8)), Expr(:kw, :max_iter, 1000)), Expr(Symbol("::"), :sempp, :SEMPPExpKern), m...),quote
        
        mts = sempp.data
        isnothing(mts) && error("No data in model, can't fit")

        model = Model(Ipopt.Optimizer)
        set_optimizer_attribute(model, "tol", tol)
        set_optimizer_attribute(model, "max_iter", max_iter)

        θ = params(sempp)
        markdens = θ[:markdens]
        impact_func = θ[:impact_function]
        
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
=#