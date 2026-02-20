mutable struct VAR{T<:Real,M<:AbstractMatrix{T},V<:AbstractMatrix{T}} <: EmpiricalModel
    Yfull::M
    p::Int64
    varnames::Vector{Symbol}

    Y::V
    X::V
    k::Int64
    t::Int64

    Phi::Union{Nothing,AbstractMatrix{T}}
    Sigma::Union{Nothing,AbstractMatrix{T}}
    E::Union{Nothing,AbstractMatrix{T}}

    F::Union{Nothing,AbstractMatrix{T}}
    Q::Union{Nothing,AbstractMatrix{T}}

    Gamma::Union{Nothing,Vector{AbstractMatrix{T}}}
    Rho::Union{Nothing,Vector{AbstractMatrix{T}}}
end

# %% Constructors

function VAR(Yfull::M, p::Int64, varnames::Vector{Symbol}) where {T<:Real,M<:AbstractMatrix{T}}
    Y, X = lagmatrix(Yfull, p)
    t, k = size(Yfull)
    V = typeof(Y)
    return VAR{T,M,V}(
        Yfull, p, varnames,          #inputs
        Y, X, k, t,                  #calculated
        nothing, nothing, nothing,   #estimate!
        nothing, nothing,            #companion!
        nothing, nothing             #other functions
    )
end

function VAR(Yfull::AbstractMatrix{T}, p::Int64) where {T<:Real}
    varnames = [Symbol("var$i") for i in 1:size(Yfull, 2)]
    return VAR(Yfull, p, varnames)
end

function VAR(Ytable, p::Int)
    Tables.istable(Ytable) ||
        throw(ArgumentError("Yfull must satisfy the Tables.jl interface"))

    Yfull = Ytable |> Tables.matrix
    names = Ytable |> Tables.columnnames .|> Symbol
    return VAR(Yfull, p, names)
end

function VAR(y::Vector, p::Int)
    return VAR(reshape(y, :, 1), p, [:y])
end

Phi(v::VAR) = v.Phi
Sigma(v::VAR) = v.Sigma
E(v::VAR) = v.E
F(v::VAR) = v.F
Q(v::VAR) = v.Q
Gamma(v::VAR) = v.Gamma
Rho(v::VAR) = v.Rho
Y(v::VAR) = v.Y
X(v::VAR) = v.X
k(v::VAR) = v.k
p(v::VAR) = v.p
t(v::VAR) = v.t
