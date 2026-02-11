mutable struct VAR{T} <: EmpiricalModel
    Yfull::Matrix{T}
    p::Int
    # varnames::Vector{Symbol}

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

    Gamma::Union{Nothing,Vector{Matrix{T}}}
    Rho::Union{Nothing,Vector{Matrix{T}}}
end



# %% Constructors

function VAR(Yfull::AbstractMatrix{T}; p::Int, varnames::Vector{Symbol}) where {T<:Real}
    Y, X = lagmatrix(Yfull, p)
    T, k = size(Yfull)
    return VAR(
        Yfull, p, varnames,                 #inputs
        Y, X, k, T, T - p,                    #calculated
        nothing, nothing, nothing, # nothing, #estimate!
        nothing, nothing, nothing,          #companion!
        nothing, nothing                    #other functions
    )
end

function VAR(Ytable; p::Int)
    Tables.istable(Ytable) ||
        throw(ArgumentError("Yfull must satisfy the Tables.jl interface"))

    Yfull = Ytable |> Tables.matrix
    names = Ytable |> Tables.columnnames .|> Symbol
    return VAR(Yfull; p=p, varnames=names)
end

function VAR(y::Vector; p::Int)
    return VAR(reshape(y, :, 1); p=p, varnames=[:y])
end
