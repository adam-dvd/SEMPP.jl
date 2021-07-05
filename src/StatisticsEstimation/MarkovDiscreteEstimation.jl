"""
    consecutive(k::Integer, sepp::SEPP)::AbstractVector{<:Real}

Compute `p_1,...,p_k` where `p_i` is the probability of an event on `i` knowing there has been events on `1,...,i-1`.
"""
function consecutive end

function consecutive(k::Integer, sepp::SEPPExpKern)::AbstractVector{<:Real}
    
    μ = sepp.μ
    ϕ = sepp.ϕ
    γ = sepp.γ

    times = 1:k
    ts = TimeSeries(times)

    vol = volfunc(times, ts, γ)
    intens = μ .+ ϕ .* vol       
    prob = 1 .- exp.(-intens )

    return prob    
end

#= TODO consecutive qui prend en compte une fonction impact pour les marques
function consecutive(k::Int, sempp::SEMPPExpKern)::AbstractVector{<:Real}
    
    μ = sempp.μ
    ϕ = sempp.ϕ
    γ = sempp.γ
    δ = sempp.δ
    markdens = sempp.markdens
    ξ = sempp.ξ
    α = sempp.α
    β = sempp.β
    κ = sempp.κ

    times = 1:k
    ts = TimeSeries(times)

    vol = volfunc(times, ts, γ)
    intens = μ .+ ϕ .* vol       
    prob = 1 .- exp.(-intens )

    return prob    
end
=#


"""
    markov_expected_run_length(sepp::SEPP, n::Integer)::Real

Compute the expected run length for the model based on a markov chain.

See Li2021 4.2.
"""
function markov_expected_run_length(sepp::SEPP, n::Integer)::Real # seulement pour SEPPExpKern pour l'instant
    p = consecutive(n,sepp)
    return (1 + sum(cumprod(p)[1:(n-1)]) + prod(p)/(1-p[n])) / (1 - sum((1-p) .* cumprod(hcat([1,1],p[2:n-1]))))
end