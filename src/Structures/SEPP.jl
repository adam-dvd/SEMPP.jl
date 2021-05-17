"""
    SEPP

Abstract type for the self-exciting point process model.
"""
abstract type SEPP end


"""
    baseline(sepp::SEPP)::Real 

Get the baseline of an SEPP.
"""
baseline(sepp::SEPP)::Real = sepp.μ


"""
    exp_decay(sepp::SEPP)::Real 

Get the exponential decay of an SEPP.

For now all SEPP models in that package have an exponential kernel, hence a method for any SEPP.
"""
exp_decay(sepp::SEPP)::Real = sepp.γ


"""
    self_excitement_factor(sepp::SEPP)::Real 

Get the factor by which the volatility function gets multiplied before it is added to the baseline in the rate.
"""
self_excitement_factor(sepp::SEPP)::Real = sepp.ϕ


"""
    params(sepp::SEPP)::AbstractDict

Get the paramaters of the model.
"""
function params end


include("SEPP/SEPPExpKern.jl")
include("SEPP/SEMPPExpKern.jl")