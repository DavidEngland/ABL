#!/usr/bin/env julia

# SHEBA-Optimized Ultraspherical Run
# Specifically tuned for Highly Stable Nocturnal Boundary Layers (HSNBL)
# Usage: julia sheba_ultra.jl data/sheba_flux.csv output/sheba_results

using CSV, DataFrames, LinearAlgebra, Statistics, LsqFit

# Optional plotting for diagnostics; script still runs without it.
const HAVE_MAKIE = try
    @eval using CairoMakie
    true
catch
    false
end

# ----------------------------- Settings ----------------------------------------

# Set true for strict MOST neutral consistency at zeta -> 0+.
const FIX_NEUTRAL_LIMIT = true

# Weighted fitting helps avoid overfitting noisy high-zeta tails.
const USE_STABILITY_WEIGHTS = true
const WEIGHT_ZETA_REF = 1.0
const WEIGHT_MIN = 0.10

# Hyperparameter grid for low-cost tuning.
const ALPHA_XI_GRID = [0.3, 0.5, 0.8]
const LAMBDA_STAR_GRID = [0.5, 0.75, 1.0]
const RIDGE_GRID = [1e-4, 1e-3, 1e-2, 5e-2]
const NMAX_GRID = [3, 4, 5, 6]

# Blocked CV in zeta space to avoid leakage across neighboring stability states.
const KFOLD = 5

# ----------------------------- Models ------------------------------------------

"""
Grachev et al. (2007) baseline for SHEBA data.
Optimized for stable conditions where traditional MOST fails.
"""
function phi_sheba_baseline(zeta, p)
    # p[1]: Neutral limit (≈1.0)
    # p[2]: 'a' constant (≈5.0)
    # p[3]: 'b' constant (≈1.1)
    # Formula: 1 + a*zeta*(1+zeta)^(1/3) / (1 + b*zeta)
    return p[1] .+ (p[2] .* zeta .* (1.0 .+ zeta).^(1/3)) ./ (1.0 .+ p[3] .* zeta)
end

"""
Grachev-style baseline with fixed neutral limit (p1 = 1.0).
This enforces strict MOST consistency at neutral conditions.
"""
function phi_sheba_baseline_fixedneutral(zeta, p)
    return 1.0 .+ (p[1] .* zeta .* (1.0 .+ zeta).^(1/3)) ./ (1.0 .+ p[2] .* zeta)
end

"""
Log-mapping for HSNBL.
Compresses the high-stability range (zeta up to 100) into the [-1, 1] interval
so Gegenbauer polynomials can resolve features across decades of stability.
"""
function xi_map_sheba(zeta, alpha_xi)
    return tanh.(alpha_xi .* log1p.(zeta))
end

function gegenbauer_eval(n::Int, lambda_star::Float64, x::Vector{Float64})
    if n == 0 return ones(length(x)) end
    if n == 1 return 2.0 * lambda_star .* x end

    c_nm1 = ones(length(x))
    c_n = 2.0 * lambda_star .* x
    for k in 1:(n - 1)
        c_np1 = (2.0 * (k + lambda_star) .* x .* c_n .- (k + 2.0 * lambda_star - 1.0) .* c_nm1) ./ (k + 1.0)
        c_nm1, c_n = c_n, c_np1
    end
    return c_n
end

function gegenbauer_design(xi::Vector{Float64}, lambda_star::Float64, nmax::Int)
    A = Matrix{Float64}(undef, length(xi), nmax + 1)
    for n in 0:nmax
        A[:, n + 1] = gegenbauer_eval(n, lambda_star, xi)
    end
    return A
end

rmse(y, yhat) = sqrt(mean((y .- yhat) .^ 2))

"""
Downweights extreme stable points while preserving a floor weight.
"""
function stability_weights(zeta::Vector{Float64}; zeta_ref::Float64=1.0, w_min::Float64=0.1)
    w = 1.0 ./ (1.0 .+ (zeta ./ zeta_ref).^2)
    return max.(w, w_min)
end

"""
Weighted ridge solve for residual coefficients:
  c = arg min ||W^(1/2)(A c - r)||^2 + alpha_reg ||c||^2
"""
function stable_solve_weighted(A::Matrix{Float64}, res::Vector{Float64}, w::Vector{Float64}, alpha_reg::Float64)
    W = Diagonal(w)
    ncoef = size(A, 2)
    lhs = transpose(A) * W * A + alpha_reg * I(ncoef)
    rhs = transpose(A) * W * res
    return lhs \ rhs
end

function stable_solve_unweighted(A::Matrix{Float64}, res::Vector{Float64}, alpha_reg::Float64)
    ncoef = size(A, 2)
    return (transpose(A) * A + alpha_reg * I(ncoef)) \ (transpose(A) * res)
end

"""
Create blocked fold indices after sorting by zeta.
"""
function make_blocked_folds(zeta::Vector{Float64}, k::Int)
    n = length(zeta)
    p = sortperm(zeta)
    fold_sizes = fill(div(n, k), k)
    for i in 1:rem(n, k)
        fold_sizes[i] += 1
    end
    folds = Vector{Vector{Int}}(undef, k)
    pos = 1
    for i in 1:k
        stop = pos + fold_sizes[i] - 1
        folds[i] = p[pos:stop]
        pos = stop + 1
    end
    return folds
end

"""
K-fold CV score for ultraspherical residual fit.
"""
function cv_score_ultra(zeta::Vector{Float64}, res::Vector{Float64}, alpha_xi::Float64, lambda_star::Float64, nmax::Int, alpha_reg::Float64)
    k = min(KFOLD, length(zeta) - 1)
    if k < 2
        error("Need at least 3 samples for blocked CV scoring.")
    end
    folds = make_blocked_folds(zeta, k)
    scores = Float64[]
    for i in eachindex(folds)
        val_idx = folds[i]
        train_parts = vcat(folds[1:(i - 1)], folds[(i + 1):end])
        train_idx = reduce(vcat, train_parts)

        z_tr = zeta[train_idx]
        z_va = zeta[val_idx]
        r_tr = res[train_idx]
        r_va = res[val_idx]

        A_tr = gegenbauer_design(xi_map_sheba(z_tr, alpha_xi), lambda_star, nmax)
        A_va = gegenbauer_design(xi_map_sheba(z_va, alpha_xi), lambda_star, nmax)

        c = if USE_STABILITY_WEIGHTS
            w_tr = stability_weights(z_tr; zeta_ref=WEIGHT_ZETA_REF, w_min=WEIGHT_MIN)
            stable_solve_weighted(A_tr, r_tr, w_tr, alpha_reg)
        else
            stable_solve_unweighted(A_tr, r_tr, alpha_reg)
        end

        rhat_va = A_va * c
        push!(scores, rmse(r_va, rhat_va))
    end
    return mean(scores)
end

# ----------------------------- Runner ------------------------------------------

function main()
    if length(ARGS) < 2
        println("Usage: julia sheba_ultra.jl <input_csv> <output_prefix>")
        return
    end

    input_csv, out_prefix = ARGS[1], ARGS[2]
    df = CSV.read(input_csv, DataFrame)

    # Filter for Stable Conditions (zeta > 0)
    df = filter(row -> !isnan(row.zeta) && row.zeta > 0 && !isnan(row.phi_obs), df)
    zeta = Vector{Float64}(df.zeta)
    y = Vector{Float64}(df.phi_obs)

    if length(y) < max(20, KFOLD + 1)
        error("Need at least $(max(20, KFOLD + 1)) stable samples after filtering for robust fitting/CV.")
    end

    # 1. Fit SHEBA baseline
    if FIX_NEUTRAL_LIMIT
        p0 = [5.0, 1.1]  # [a, b]
        fit = curve_fit(phi_sheba_baseline_fixedneutral, zeta, y, p0)
        p_base = fit.param
        yhat_base = phi_sheba_baseline_fixedneutral(zeta, p_base)
    else
        p0 = [1.0, 5.0, 1.1]  # [neutral, a, b]
        fit = curve_fit(phi_sheba_baseline, zeta, y, p0)
        p_base = fit.param
        yhat_base = phi_sheba_baseline(zeta, p_base)
    end

    # 2. Residual fit target
    res = y .- yhat_base

    # 3. Hyperparameter search (blocked CV on residual prediction)
    best = nothing
    best_cv = Inf
    for alpha_xi in ALPHA_XI_GRID
        for lambda_star in LAMBDA_STAR_GRID
            for alpha_reg in RIDGE_GRID
                for nmax in NMAX_GRID
                    score = cv_score_ultra(zeta, res, alpha_xi, lambda_star, nmax, alpha_reg)
                    if score < best_cv
                        best_cv = score
                        best = (alpha_xi=alpha_xi, lambda_star=lambda_star, alpha_reg=alpha_reg, nmax=nmax)
                    end
                end
            end
        end
    end

    # 4. Final fit with best hyperparameters
    xi = xi_map_sheba(zeta, best.alpha_xi)
    A = gegenbauer_design(xi, best.lambda_star, best.nmax)
    coeffs = if USE_STABILITY_WEIGHTS
        w = stability_weights(zeta; zeta_ref=WEIGHT_ZETA_REF, w_min=WEIGHT_MIN)
        stable_solve_weighted(A, res, w, best.alpha_reg)
    else
        stable_solve_unweighted(A, res, best.alpha_reg)
    end
    yhat_ultra = yhat_base .+ A * coeffs

    # 5. Output metrics
    println("--- SHEBA HSNBL Fit Results ---")
    println("Baseline (Grachev) RMSE: ", rmse(y, yhat_base))
    println("Ultra-Corrected RMSE:    ", rmse(y, yhat_ultra))
    println("Ultra residual CV RMSE:  ", best_cv)
    if FIX_NEUTRAL_LIMIT
        println("Base Params [a, b] with neutral fixed at 1.0: ", p_base)
    else
        println("Base Params [neutral, a, b]: ", p_base)
    end
    println("Best hyperparameters: ", best)

    # 6. Diagnostics for heteroscedastic tail behavior
    res_base = y .- yhat_base
    res_ultra = y .- yhat_ultra
    println("Residual std (baseline): ", std(res_base))
    println("Residual std (ultra):    ", std(res_ultra))

    # Save results
    res_df = DataFrame(
        zeta=zeta,
        obs=y,
        baseline=yhat_base,
        ultra=yhat_ultra,
        residual_baseline=res_base,
        residual_ultra=res_ultra,
    )
    CSV.write("$(out_prefix)_sheba_pred.csv", res_df)

    par_df = DataFrame(
        fixed_neutral=[FIX_NEUTRAL_LIMIT],
        baseline_param_1=[p_base[1]],
        baseline_param_2=[p_base[2]],
        baseline_param_3=[FIX_NEUTRAL_LIMIT ? NaN : p_base[3]],
        alpha_xi=[best.alpha_xi],
        lambda_star=[best.lambda_star],
        alpha_reg=[best.alpha_reg],
        nmax=[best.nmax],
        use_weights=[USE_STABILITY_WEIGHTS],
        weight_zeta_ref=[WEIGHT_ZETA_REF],
        weight_min=[WEIGHT_MIN],
        cv_rmse=[best_cv],
    )
    CSV.write("$(out_prefix)_sheba_params.csv", par_df)

    # Optional figures for quick quality check.
    if HAVE_MAKIE
        fig1 = Figure(resolution=(900, 540))
        ax1 = Axis(fig1[1, 1], xlabel="zeta", ylabel="phi", title="SHEBA Baseline vs Ultraspherical")
        scatter!(ax1, zeta, y, markersize=7, color=(:black, 0.55), label="obs")
        scatter!(ax1, zeta, yhat_base, markersize=5, color=(:orangered, 0.45), label="baseline")
        scatter!(ax1, zeta, yhat_ultra, markersize=5, color=(:seagreen, 0.45), label="ultra")
        axislegend(ax1, position=:lt)
        save("$(out_prefix)_sheba_fit.png", fig1)

        fig2 = Figure(resolution=(900, 540))
        ax2 = Axis(fig2[1, 1], xlabel="zeta", ylabel="residual", title="Residuals vs zeta")
        scatter!(ax2, zeta, res_base, markersize=6, color=(:orangered, 0.5), label="baseline residual")
        scatter!(ax2, zeta, res_ultra, markersize=6, color=(:seagreen, 0.5), label="ultra residual")
        hlines!(ax2, [0.0], color=:black, linewidth=1.2)
        axislegend(ax2, position=:rt)
        save("$(out_prefix)_sheba_residuals_vs_zeta.png", fig2)
    end

    println("\nResults saved to $(out_prefix)_sheba_pred.csv")
    println("Parameters saved to $(out_prefix)_sheba_params.csv")
    if HAVE_MAKIE
        println("Saved figures: $(out_prefix)_sheba_fit.png, $(out_prefix)_sheba_residuals_vs_zeta.png")
    else
        println("Plotting skipped: CairoMakie not installed.")
    end
end

main()
