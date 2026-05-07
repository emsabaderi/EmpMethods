# %% module def, imports

# import ShiftedArrays as SA

# using Random

# import ComponentArrays as CA

# %% helper functions

function lagmatrix(Y::AbstractMatrix{T}, p::Int) where {T<:Real}
    t, k = size(Y)
    @assert t > p "More lags than observations"
    n = t - p
    X = Matrix{T}(undef, n, k * p)
    for i in 1:p
        X[:, (1:k).+k*(i-1)] .= Y[p-i+1:t-i, :]
    end
    return Y[p+1:end, :], X
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
    T = eltype(var.Y)

    lagnames = Vector{Symbol}(undef, p)
    Phis = Vector{AbstractMatrix{T}}(undef, p)
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
        F[1:k, (1:k).+i*k] = (@view Phi[(1:k).+i*k, :])'
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
    T = eltype(var.Y)

    Γ0 = lyapd(F, Q)

    # Preallocate
    Gamma = Vector{Matrix{T}}(undef, H + 1)
    Gamma[1] = Γ0[1:k, 1:k]

    Fh = Matrix{T}(I, kp, kp)
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

function autocov!(svar::SVAR; H=20)
    autocov!(svar.var; H=H)
end

# %%

function autocorr!(var::VAR)
    Γ0 = var.Gamma[1]
    stds = sqrt.(diag(Γ0))
    Dinv = Diagonal(1 ./ stds)

    var.Rho = [Dinv * Γh * Dinv for Γh in var.Gamma]
    return nothing
end

function autocorr!(svar::SVAR)
    autocorr!(svar.var)
end

# %%

function shortrun!(svar::SVAR)
    Sigma = svar.var.Sigma
    E = svar.var.E
    svar.A = cholesky(Sigma).L
    svar.U = (svar.A \ E')'
    return nothing
end

function longrun!(svar::SVAR)
    all(abs.(eigvals(svar.var.F)) .< 1) ||
        throw(ArgumentError("VAR object is not stationary. Long-run identification requires eigenvalues of VAR companion matrix to be within the unit circle"))
    Σ = Sigma(svar)
    blocks = phiblocks(svar.var)
    PhiSum = sum(blocks[lag] for lag in keys(blocks))
    IminPhi = I - PhiSum
    Q = IminPhi \ (Σ / IminPhi')
    Q = (Q + Q') / 2
    svar.A = IminPhi * cholesky(Q).L
    svar.U = (svar.A \ E(svar)')'
    return nothing
end

# %%

function imresp(v::VAR; horizon::Int=20)
    Φ = v.Phi
    Fmat = v.F
    k = v.k
    kp = size(Fmat, 1)
    T = eltype(Φ)

    IRFs = Array{T}(undef, k, k, horizon + 1)
    Fh = Matrix{T}(I, kp, kp)
    for h in 0:horizon
        IRFs[:, :, h+1] = Fh[1:k, 1:k]
        Fh = Fh * Fmat
    end

    return IRFs
end

function imresp(s::SVAR; horizon::Int=20)
    isnothing(s.A) && error("SVAR not identified. Run shortrun! or longrun! first.")
    IRFs = imresp(s.var; horizon=horizon)
    for h in 1:size(IRFs, 3)
        IRFs[:, :, h] = IRFs[:, :, h] * s.A
    end
    return IRFs
end
