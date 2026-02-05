# EmpMethods

[![Build Status](https://github.com/Emad/EmpMethods.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/Emad/EmpMethods.jl/actions/workflows/CI.yml?query=branch%3Amaster)

EmpMethods/
  src/
    EmpMethods.jl
    var/
      VAR.jl
      estimation.jl
      irf.jl
      utils.jl
    svar/
      SVAR.jl
      identification.jl
      irf.jl
    dsge/
      DSGE.jl
      solution.jl
      simulation.jl
    bayes/
      BayesianVAR.jl
      priors.jl
      gibbs.jl
      
You want every model type to follow the same pattern:
A struct that stores:
- model parameters
- estimated objects
- data
- metadata
A constructor that performs estimation or setup
Methods that operate on the struct
- irf(model, horizon)
- forecast(model, steps)
- simulate(model, ...)
- loglik(model)
- plot(model)
This gives you a uniform API across all models.
