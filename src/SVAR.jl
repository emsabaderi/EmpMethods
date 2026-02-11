mutable struct SVAR{T} <: EmpiricalModel
    Yfull::Matrix{T}
    p::Int

    Y::Matrix{T}
    X::Matrix{T}
    k::Int
    T::Int
    Te::Int

    Phi::Union{Nothing,Matrix{T}}
    Sigma::Union{Nothing,Matrix{T}}
    # coeftable::Union{Nothing,DataFrames.DataFrame}
    E::Union{Nothing,Matrix{T}}

    Phi_blocks::Union{Nothing,Vector{Matrix{T}}}
    F::Union{Nothing,Matrix{T}}
    Q::Union{Nothing,Matrix{T}}

    A::Union{Nothing,Matrix{T}}
    C1::Union{Nothing,Matrix{T}}
    scheme::Symbol
end

function SVAR(Yfull::AbstractMatrix{T}, A::AbstractMatrix{T}; p::Int, scheme::Symbol=:none) where {T<:Real}
    Y, X = lagmatrix(Yfull, p)
    T, k = size(Yfull)


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
