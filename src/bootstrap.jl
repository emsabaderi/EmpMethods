function bootstrap(var::VAR; N::Int=1000, horizon::Int=20,
    scheme::Symbol=:short_run, level::Float64=0.90)

    k = var.k
    p = var.p
    Tobs = size(var.Y, 1)

    # 1. Original estimates
    Phi_hat = var.Phi              # (k*p)×k
    eps_hat = var.residuals        # T×k residuals
    y0 = var.Y[1:p, :]             # initial conditions (p×k)

    # Storage: (k variables, k shocks, horizon+1, N draws)
    IRF_store = zeros(eltype(var.Y), k, k, horizon + 1, N)

    # Helper: simulate VAR given Phi, shocks, and initial conditions
    function simulate_var(Phi::Matrix{T}, shocks::Matrix{T},
        y0::Matrix{T}, k::Int, p::Int) where {T}

        Tsim = size(shocks, 1)
        Ysim = zeros(Tsim, k)

        # fill initial conditions
        for i in 1:p
            Ysim[i, :] .= y0[i, :]
        end

        # simulate forward
        for t in p+1:Tsim
            x = vec(reverse(Ysim[t-p:t-1, :], dims=1))  # length k*p
            Ysim[t, :] = Phi' * x + shocks[t, :]       # k×(k*p) * (k*p) = k
        end

        return Ysim
    end

    # 2. Bootstrap loop
    for n in 1:N
        # 2.1 draw residuals with replacement
        idx = rand(1:Tobs, Tobs)
        eps_boot = eps_hat[idx, :]

        # 2.2 simulate new data using Phi_hat, eps_boot, and y0
        Y_boot = simulate_var(Phi_hat, eps_boot, y0, k, p)

        # 2.3 re-estimate VAR on Y_boot
        var_boot = VAR(Y_boot; p=p, varnames=[:x1, :x2])
        estimate!(var_boot)
        companion!(var_boot)

        # 2.4 identify SVAR, compute IRFs, store
        svar_boot = SVAR(var_boot; scheme=scheme)
        IRF_store[:, :, :, n] = irf(svar_boot; horizon=horizon)
    end

    # 3. Percentiles across bootstrap draws
    α = (1 - level) / 2
    lower = similar(IRF_store, k, k, horizon + 1)
    median = similar(IRF_store, k, k, horizon + 1)
    upper = similar(IRF_store, k, k, horizon + 1)

    for i in 1:k, j in 1:k, h in 1:(horizon+1)
        v = view(IRF_store, i, j, h, :)
        lower[i, j, h] = Statistics.quantile(v, α)
        median[i, j, h] = Statistics.quantile(v, 0.5)
        upper[i, j, h] = Statistics.quantile(v, 1 - α)
    end

    return lower, median, upper
end
