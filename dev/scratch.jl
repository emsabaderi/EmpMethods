# %% Development scratch file for testing EmpMethods functions
# Workflow: Edit functions in src/, changes auto-reload via Revise, test here

using Revise
using BenchmarkTools
using Random
using DataFrames
using Statistics


# %%

using EmpMethods

# ============================================================================
# %% Test controlled inputs
# ============================================================================

Random.seed!(654)
p, k = 3, 2;
Y = rand(Float64, 100, k);

# %%

@btime lagmatrix1(Y, p)

# %%

@btime lagmatrix2(Y, p)

# %%
myvar = VAR(Y; p=p, varnames=[:var1, :var2]);
# %%
@btime estimate!(myvar)

# %%

myvar.Phi

# %%
# ============================================================================
# Test random inputs
# ============================================================================

# Random.seed!(42)
# # Your random test cases here

# ============================================================================
# Debugging / Development
# ============================================================================

# Place breakpoints or detailed testing here

# %%
