# %% Development scratch file for testing EmpMethods functions
# Workflow: Edit functions in src/, changes auto-reload via Revise, test here

using Revise
using EmpMethods
using BenchmarkTools
using Random
using DataFrames
using Statistics

# ============================================================================
# %% Test controlled inputs
# ============================================================================

Random.seed!(654)
p,k = 3,2 ;
Y = rand(Float64, 100, k) ;
tvar = VAR(Y; p=p, varnames=[:dy, :inf]) ;
estimate!(tvar) ;
companion!(tvar) ;

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
