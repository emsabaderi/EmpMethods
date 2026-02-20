# %%

module EmpMethods

using Statistics, DataFrames, ShiftedArrays, StaticArrays, LinearAlgebra, ComponentArrays, MatrixEquations

abstract type EmpiricalModel end

include("VAR.jl")
include("SVAR.jl")
include("utils.jl")

export VAR, Phi, Sigma, E, F, Q, Gamma, Rho, Y, X, k, p, t
export SVAR
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
