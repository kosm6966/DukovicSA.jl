module DukovicSA

import DataFrames
import CSV
using Optim
using DelimitedFiles
using Interpolations
using Random
using Parameters
using SpecialFunctions
using StaticArrays

include("parameters.jl")
include("data_functions.jl")
include("rfKMC_HT.jl")
include("rfKMC_RC.jl")
include("rfKMC_INIT.jl")
include("optimize_function.jl")
# include("optim/api.jl")
# include("optim/interface.jl")
# include("optim/assess_convergence.jl")
# include("optim/types.jl")
# include("optim/samin.jl")
# include("optim/update.jl")
# include("optim/maximize.jl")
# include("optim/generic.jl")


export StFunction
export Data
export KMC_HT
export KMC_RC
export run_optim
export saveascsv
export makerandp0file
export fitKMC

end
