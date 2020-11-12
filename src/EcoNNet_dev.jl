"""
General convention on ordering variables:
	1) Order in chronological order (i.e. lagged, current, lead)
	2) Order alphabetically within each period
"""

# Import dependencies
using Distributions, Parameters, LinearAlgebra, QuantEcon
using Flux, Distributed, SharedArrays, JLD, NLsolve, ProgressMeter
using Flux: throttle, params, mse, glorot_uniform, @epochs
using IterTools: ncycle, NCycle
using NLsolve: SolverResults
using DataFrames, Random, Atom, Juno, DataFramesMeta
using Plots, Plots.PlotMeasures



include("structs.jl")
include("object_management.jl")
include("stochastic.jl")
include("aux_fns.jl")
include("beliefs.jl")
include("state_df.jl")
include("step_fns.jl")
include("grid_approach.jl")
include("sim_approach.jl")
include("plotting.jl")
