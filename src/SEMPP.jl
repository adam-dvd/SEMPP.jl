module SEMPP

using Distributions
using EGPD
using JuMP
using Ipopt

include("Types.jl")
include("DiscreteHawkesProcess.jl")

export 
    # data types
    PP, PointProcess, MarkedPointProcess

    # data types methods
    end_time, start_time, ground_process

    # model types
    SEPP, DiscreteSEPPExpKern, DiscreteSEMPPExpKern

    # model types methods
    baseline, params, exp_decay, self_excitement_factor, params 
    lin_coeff_impact, marks_scale_params, shape, decay

    # model manipulation
    negloglik, fit!

end
