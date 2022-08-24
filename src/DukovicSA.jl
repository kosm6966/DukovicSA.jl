module DukovicSA

import DataFrames
import CSV
using Optim
using DelimitedFiles
using Interpolations
using Distributions
using StatsBase
using Random
using Parameters
using StaticArrays
using SpecialFunctions
using StatsPlots
using Plots

include("parameters.jl")
include("data_functions.jl")
include("rfKMC_HT.jl")
include("rfKMC_RC.jl")
include("rfKMC_INIT.jl")
include("optimize_function.jl")

export StFunction
export Data
export KMC_HT
export KMC_RC
export run_optim
export saveascsv
export makerandp0file
export fitKMC
export StFunction
export randomdistribution

end
