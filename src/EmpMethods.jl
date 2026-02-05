# %%

module EmpMethods

import Statistics
import DataFrames
import LinearAlgebra: I, Diagonal, kron, sqrt, diag, cholesky, diagm
import ShiftedArrays
import ComponentArrays
import Tables

abstract type EmpiricalModel end

include("VAR.jl")
include("SVAR.jl")
include("utils.jl")

export VAR
export SVAR
export lagmatrix
export estimate!
export companion!
export autocov!
export autocorr!
export longrun_impact
export irf
export bootstrap

end

# %%
