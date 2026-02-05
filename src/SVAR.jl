mutable struct SVAR{T} <: EmpiricalModel
    var::VAR{T}
    A::Union{Nothing,Matrix{T}}
    C1::Union{Nothing,Matrix{T}}
    scheme::Symbol
end

function SVAR(var::VAR{T}; scheme::Symbol=:none) where {T}
    C1 = nothing
    if scheme == :none
        A = nothing
    elseif scheme == :short_run
        A = cholesky(var.Sigma).L
    elseif scheme == :long_run
        C1 = longrun_impact(var)
        A = cholesky(C1 * var.Sigma * C1').L
    else
        error("Unknown SVAR identification scheme: $scheme")
    end

    return SVAR{T}(var, A, C1, scheme)
end

function SVAR(var::VAR{T}, A::Matrix{T}; scheme::Symbol=:user) where {T}
    return SVAR{T}(var, A, nothing, scheme)
end
