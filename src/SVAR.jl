mutable struct SVAR{T<:Real,AType<:AbstractMatrix{T}} <: EmpiricalModel
    var::VAR{T}
    A::AType
    impact::Union{Nothing,Matrix{T}}  # structural impact matrix (A^{-1}, A^{-1}P, etc.)
end

function SVAR(var::VAR{T}, A::AbstractMatrix) where {T}
    Afloat = Array{T}(A)
    SVAR{T,typeof(Afloat)}(var, Afloat, nothing)
end

function SVAR(Yfull::AbstractMatrix{T}, p::Int, names, A::AbstractMatrix) where {T}
    var = VAR(Yfull, p, names)
    estimate!(var)
    companion!(var)
    autocov!(var)
    autocorr!(var)
    SVAR(var, A)
end

function SVAR(Ytable, p::Int, A::AbstractMatrix) where {T}
    Tables.istable(Ytable) ||
        throw(ArgumentError("Yfull must satisfy the Tables.jl interface"))

    Yfull = Tables.matrix(Ytable)
    names = Tables.columnnames(Ytable) .|> Symbol
    SVAR(Yfull, p, names, A)
end

Phi(s::SVAR) = s.var.Phi
Sigma(s::SVAR) = s.var.Sigma
E(s::SVAR) = s.var.E
F(s::SVAR) = s.var.F
Q(s::SVAR) = s.var.Q
Y(s::SVAR) = s.var.Y
X(s::SVAR) = s.var.X
k(s::SVAR) = s.var.k
p(s::SVAR) = s.var.p
t(s::SVAR) = s.var.t
