struct Bootstrap
    var::Union{VAR, SVAR}
    horizon::Int
    seed::Int
end

function Bootstrap(svar::SVAR, horizon::Int, seed::Int)
    Random.seed!(seed)

    var = svar.var
    estimate!(var)
    Phi0 = var.Phi
    residuals0 = var.residuals
