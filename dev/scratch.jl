# %% Development scratch file for testing EmpMethods functions
# Workflow: Edit functions in src/, changes auto-reload via Revise, test here

using Revise
using BenchmarkTools
using Random
using DataFrames
using Statistics
using StaticArrays

# %%

using EmpMethods

# ============================================================================
# %% Test controlled inputs
# ============================================================================

Random.seed!(634)
p, k = 3, 2;
Y = rand(-1:0.01:1, 100, k);
# %%
myvar = VAR(Y, p, [:var1, :var2]);
estimate!(myvar)
companion!(myvar)
autocov!(myvar)
autocorr!(myvar)
# %%
@btime myvar = VAR(Y, p, [:var1, :var2]);
# %%
@btime estimate!(myvar);
# %%
@btime companion!(myvar);
# %%
@btime phiblocks(myvar);
# %%
@btime autocov!(myvar);
# %%
@btime autocorr!(myvar);
# %%
mysvar = SVAR(Y, p, [:var1, :var2], [1 0; 1 1])
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
