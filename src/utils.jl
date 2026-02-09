# %% module def, imports

# import ShiftedArrays as SA

# using Random

# import ComponentArrays as CA

# %% helper functions

"""
function lagmatrix(Y::Matrix, p::Int)
    ...
    return X, Ytrim
end

Given a matrix of time series variables, outputs the p - 1 scaled dependent variable matrix and the progressively-lagged (from t - 1 to t - p) explanatory variable matrix as a tuple.

Ensure data is fed in ascending order (t = 1 is earliest date).
"""
function lagmatrix(Y::Matrix{Float64}, p::Int)
    Ymat = ndims(Y) == 1 ? reshape(Y, :, 1) : Y
    T, k = size(Ymat)
    @assert T > p "More lags than observations"

    Ytrim = Ymat[p+1:end, :]

    X = Matrix{eltype(Y)}(hcat([ShiftedArrays.lag(Y, lag)[p+1:T, :] for lag in 1:p]...))

    return (Ytrim=Ytrim, X=X)
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
    T = var.T
    lower = hcat(Matrix(I(k * p - k)), zeros(Float64, (k * p - k, k * p - k)))
    names = var.varnames
    phinames = Symbol[]
    for lag in 1:var.p
        for nm in names
            push!(phinames, Symbol("$(nm)_L$(lag)"))
        end
    end
    Phi = X \ Y

    coeftable = hcat(
        DataFrames.DataFrame(term=phinames),
        DataFrames.DataFrame(Phi, names)
    )

    R = Y .- X * Phi
    Sigma = R'R ./ (T - p * k)

    # F = vcat(Phi', lower)

    var.Phi = Phi
    var.Sigma = Sigma
    var.coeftable = coeftable
    var.residuals = R
    # var.F = F

    return nothing
end

# %%
"""
function companion!(var::VAR)
    ...
    return nothing
end
Given a var element, calculates its companion-form F matrix and Q matrix.

ensure "estimate!" has been called first.
"""

function companion!(var::VAR{T}) where {T}
    k = var.k
    p = var.p
    Phi = var.Phi

    Phi_blocks = [Phi[(i-1)*k+1:i*k, :] for i in 1:p]

    F = zeros(T, k * p, k * p)

    F[1:k, 1:k*p] = hcat(Phi_blocks...)

    if p > 1
        F[k+1:end, 1:k*(p-1)] = I(k * (p - 1))
    end

    Q = zeros(T, k * p, k * p)
    Q[1:k, 1:k] .= var.Sigma

    var.Phi_blocks = Phi_blocks
    var.F = F
    var.Q = Q

    return nothing
end


# %%

function autocov!(var::VAR; H=20)
    F = var.F
    Q = var.Q
    k = var.k
    p = var.p
    kp = k * p

    A = I(kp^2) - kron(F, F)
    Γ0_vec = A \ vec(Q)
    Γ0 = reshape(Γ0_vec, kp, kp)

    var.Gamma = Vector{Matrix{Float64}}(undef, H + 1)
    var.Gamma[1] = Γ0[1:k, 1:k]

    for h in 1:H
        Γx = F^h * Γ0
        var.Gamma[h+1] = Γx[1:k, 1:k]
    end

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

function longrun_impact(var::VAR{T}) where {T}
    Φs = var.Phi_blocks
    isnothing(Φs) && error("Phi_blocks is not computed. Run companion! on the VAR first.")

    k = size(Φs[1], 1)

    S = zero(Φs[1])
    for Φ in Φs
        S .+= Φ
    end

    C1 = inv(I(k) - S)
    return C1
end


# %%

function irf(s::SVAR{T}; horizon::Int=20) where {T}
    isnothing(s.A) && error("SVAR is not identified. A matrix is missing.")

    F = s.var.F
    A = s.A
    k = s.var.k

    IRFs = Array{T}(undef, k, k, horizon + 1)

    IRFs[:, :, 1] = A

    Fpow = I(size(F, 1))
    for h in 1:horizon
        Fpow = Fpow * F
        IRFs[:, :, h+1] = Fpow[1:k, 1:k] * A
    end

    return IRFs
end
