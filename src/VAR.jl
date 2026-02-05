mutable struct VAR{T} <: EmpiricalModel
    # inputs
    Yfull::Matrix{T}
    p::Int
    varnames::Vector{Symbol}

    # calculated
    Y::Matrix{T}
    X::Matrix{T}
    k::Int
    T::Int
    Te::Int

    # via estimate!
    Phi::Union{Nothing,Matrix{T}}
    Sigma::Union{Nothing,Matrix{T}}
    coeftable::Union{Nothing,DataFrames.DataFrame}
    residuals::Union{Nothing,Matrix{T}}

    # via companion!
    Phi_blocks::Union{Nothing,Vector{Matrix{T}}}
    F::Union{Nothing,Matrix{T}}
    Q::Union{Nothing,Matrix{T}}

    # via other functions
    Gamma::Union{Nothing,Vector{Matrix{T}}}
    Rho::Union{Nothing,Vector{Matrix{T}}}
end



# %% Constructors

function VAR(Yfull::Matrix{Float64}; p::Int, varnames::Vector{Symbol})
    Y, X = lagmatrix(Yfull, p)
    T, k = size(Yfull)
    return VAR(
        Yfull, p, varnames,                 #inputs
        Y, X, k, T, T - p,                    #calculated
        nothing, nothing, nothing, nothing, #estimate!
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
