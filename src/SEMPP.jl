module SEMPP

using Distributions
using EGPD
using JuMP
using Ipopt
using DataFrames, Gadfly
using Dates
using Extremes
using Statistics
using ForwardDiff

include("Structures.jl")
include("ParametersEstimation.jl")
include(joinpath("Plots", "PlotsData.jl"))
include(joinpath("Plots", "ModelPlots.jl"))
include(joinpath("Plots", "ValidationPlots.jl"))
include("Simulation.jl")
include("StatisticsEstimation.jl")

export 
    # data types
    TS, TimeSeries, MarkedTimeSeries,

    # data types methods
    end_time, start_time, ground_process, copy_above,

    # model types
    SEPP, SEPPExpKern, SEMPPExpKern,

    # model types methods
    baseline, params, exp_decay, self_excitement_factor, params, 
    lin_coeff_impact, marks_scale_params, shape, decay,

    # model manipulation
    discrete_negloglik, negloglik, discrete_fit!, fit!,

    # plots
    rate_plot, pp_prob_plot, marks_prob_plot, marks_qq_plot,

    # simulations
    simulation, discrete_simulation,

    # statistics estimation
    ## Markov
    markov_expected_run_length,

    ## Monte Carlo 
    monte_carlo_return_period, rQy,

    ## tail estimation
    discrete_tail_estimation, tail_estimation

end # module
