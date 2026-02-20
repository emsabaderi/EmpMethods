# %% module def, imports

# import ShiftedArrays as SA

# using Random

# import ComponentArrays as CA

# %% helper functions

function lagmatrix(Y::Matrix{Float64}, p::Int)
    Y = ndims(Y) == 1 ? reshape(Y, :, 1) : Y
    t, k = size(Y)
    @assert t > p "More lags than observations"

    X = fill(NaN, t, k * p)

    for i in 1:p
        @views X[(i+1):t, (1:k).+k*(i-1)] .= Y[1:(t-i), :]
    end
    return @views Y[p+1:end, :], @views X[p+1:end, :]
end

# %%
"""
function estimate!(var::VAR)
    ...
    return nothing
end
Given a VAR element, estimates its parameters and variance-covariance matrix with OLS.
"""
function estimate!(var::VAR)
    X = var.X
    Y = var.Y
    k = var.k
    p = var.p
    t = var.t

    Phi = X \ Y
    E = Y .- X * Phi
    sdn = 1 / (t - p * k)

    var.Phi = Phi
    var.E = E
    var.Sigma = E'E .* sdn

    return nothing
end

# %%

function phiblocks(var::VAR)
    k = var.k
    p = var.p

    lagnames = Vector{Symbol}(undef, p)
    Phis = Vector{AbstractMatrix{Float64}}(undef, p)
    for lag in 1:p
        lagnames[lag] = Symbol("lag", lag)
        Phis[lag] = @view var.Phi[(1:k).+(lag-1)*k, :]
    end

    return ComponentArray(; zip(lagnames, Phis)...)
end

"""
function companion!(var::VAR)
    ...
    return nothing
end
Given a var element, calculates its companion-form F matrix and Q matrix.

ensure "estimate!" has been called first.
"""
function companion!(var::VAR)
    k = var.k
    p = var.p
    T = eltype(var.Y)
    kp = k * p
    Phi = var.Phi

    # F
    F = zeros(T, kp, kp)
    IF = I(kp - k)

    for i in 0:(p-1)
        F[1:k, (1:k).+i*k] = @view Phi[(1:k).+i*k, :]
    end

    if p > 1
        F[k+1:end, 1:kp-k] = IF
    end

    # Q
    Q = zeros(T, kp, kp)
    Q[1:k, 1:k] .= var.Sigma

    var.F = F
    var.Q = Q

    return nothing
end


# %%


function autocov!(var::VAR; H=20)
    F = var.F
    Q = var.Q
    k = var.k
    kp = k * var.p

    Γ0 = lyapd(F, Q)

    # Preallocate
    Gamma = Vector{Matrix{Float64}}(undef, H + 1)
    Gamma[1] = Γ0[1:k, 1:k]

    Fh = Matrix{Float64}(I, kp, kp)
    Fhtmp = similar(Fh)
    Γx = similar(Γ0)

    for h in 1:H
        mul!(Fhtmp, Fh, F)
        Fh, Fhtmp = Fhtmp, Fh
        mul!(Γx, Fh, Γ0)
        # Γx = Fh * Γ0
        Gamma[h+1] = Γx[1:k, 1:k]
    end

    var.Gamma = Gamma
    return nothing
end

# %%

function autocorr!(var::VAR)
    Γ0 = var.Gamma[1]
    stds = sqrt.(diag(Γ0))
    Dinv = Diagonal(1 ./ stds)

    var.Rho = [Dinv * Γh * Dinv for Γh in var.Gamma]
    return nothing
end

# %%

# function longrun_impact(var::VAR{T}) where {T}
#     Φs = var.Phi_blocks
#     isnothing(Φs) && error("Phi_blocks is not computed. Run companion! on the VAR first.")

#     k = size(Φs[1], 1)

#     S = zero(Φs[1])
#     for Φ in Φs
#         S .+= Φ
#     end

#     C1 = inv(I(k) - S)
#     return C1
# end


# %%

# function irf(s::SVAR{T}; horizon::Int=20) where {T}
#     isnothing(s.A) && error("SVAR is not identified. A matrix is missing.")

#     F = s.var.F
#     A = s.A
#     k = s.var.k

#     IRFs = Array{T}(undef, k, k, horizon + 1)

#     IRFs[:, :, 1] = A

#     Fpow = I(size(F, 1))
#     for h in 1:horizon
#         Fpow = Fpow * F
#         IRFs[:, :, h+1] = Fpow[1:k, 1:k] * A
#     end

#     return IRFs
# end
