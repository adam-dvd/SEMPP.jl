module SEMPP

using Distributions
using EGPD
using JuMP
using Ipopt
using DataFrames, Gadfly
using Dates

include("Structures.jl")
include("ParametersEstimation.jl")
include(joinpath("Plots", "PlotsData.jl"))
include(joinpath("Plots", "ModelPlots.jl"))
include(joinpath("Plots", "ValidationPlots.jl"))

export 
    # data types
    PP, PointProcess, MarkedPointProcess,

    # data types methods
    end_time, start_time, ground_process,

    # model types
    SEPP, SEPPExpKern, SEMPPExpKern,

    # model types methods
    baseline, params, exp_decay, self_excitement_factor, params, 
    lin_coeff_impact, marks_scale_params, shape, decay,

    # model manipulation
    discrete_negloglik, negloglik, discrete_fit!, fit!,

    # plots
    rate_plot, pp_prob_plot

end
