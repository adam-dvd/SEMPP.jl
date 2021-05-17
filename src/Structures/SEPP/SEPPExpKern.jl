"""
SEPPExpKern(μ, ϕ, γ)

SEPP with exponential kernel model type with rate λ (for all t, there is an event with probability λ(t)) viz. :
```λ(t) = μ + ϕ \\sum_{t_k < t} \\exp(γ(t-t_k))```
where `t_k` are the timestamps of all events.
"""
mutable struct SEPPExpKern <: SEPP
    data::Union{TS, Nothing}
    μ::Real
    ϕ::Real
    γ::Real


    function SEPPExpKern(data::Union{TS, Nothing} = nothing; μ::Real = rand(), ϕ::Real = rand(), γ::Real = rand()) 
        if any((μ, ϕ, γ) .< 0) 
            error("paramaters must be positive or zero")
        else 
            new(data,μ, ϕ, γ)
        end
    end

end


params(sepp::SEPPExpKern)::AbstractDict = Dict(:μ => sepp.μ, :ϕ => sepp.ϕ, :γ => sepp.γ)


