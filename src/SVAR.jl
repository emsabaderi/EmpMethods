mutable struct SVAR{T<:Real,M<:AbstractMatrix{T},V<:AbstractMatrix{T}} <: EmpiricalModel
    var::VAR{T,M,V}
    A::Union{AbstractMatrix{T},Nothing}
    U::Union{AbstractMatrix{T},Nothing}
end

function SVAR(var::VAR{T,M,V}) where {T<:Real,M<:AbstractMatrix{T},V<:AbstractMatrix{T}}
    return SVAR{T,M,V}(var, nothing, nothing)
end

function SVAR(Yfull::AbstractMatrix{T}, p::Int, names::Vector{Symbol}) where {T<:Real}
    var = VAR(Yfull, p, names)
    estimate!(var)
    companion!(var)
    autocov!(var)
    autocorr!(var)
    return SVAR(var)
end

function SVAR(Ytable, p::Int)
    Tables.istable(Ytable) ||
        throw(ArgumentError("Yfull must satisfy the Tables.jl interface"))

    Yfull = Tables.matrix(Ytable)
    names = Tables.columnnames(Ytable) .|> Symbol
    return SVAR(Yfull, p, names)
end

Phi(s::SVAR) = s.var.Phi
Sigma(s::SVAR) = s.var.Sigma
E(s::SVAR) = s.var.E
F(s::SVAR) = s.var.F
Y(s::SVAR) = s.var.Y
X(s::SVAR) = s.var.X
k(s::SVAR) = s.var.k
p(s::SVAR) = s.var.p
t(s::SVAR) = s.var.t
U(s::SVAR) = s.U
A(s::SVAR) = s.A
