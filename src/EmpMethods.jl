# %%

module EmpMethods

using Statistics, DataFrames, ShiftedArrays, StaticArrays, LinearAlgebra, ComponentArrays

abstract type EmpiricalModel end

include("VAR.jl")
# include("SVAR.jl")
include("utils.jl")

export VAR
# export SVAR
export lagmatrix
export estimate!
export phiblocks
export companion!
export autocov!
export autocorr!
# export longrun_impact
# export irf
# export bootstrap

end

# %%
