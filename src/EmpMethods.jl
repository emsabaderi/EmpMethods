# %%

module EmpMethods

using Statistics, DataFrames, ShiftedArrays, StaticArrays, LinearAlgebra, ComponentArrays, MatrixEquations

abstract type EmpiricalModel end

include("VAR.jl")
include("SVAR.jl")
include("utils.jl")

export VAR
export SVAR, Phi, Sigma, E, F, Q, Gamma, Rho, Y, X, k, p, t, A, U
export lagmatrix, estimate!, phiblocks, companion!, autocov!, autocorr!, shortrun!, longrun!
# export longrun_impact
# export irf
# export bootstrap

end

# %%
