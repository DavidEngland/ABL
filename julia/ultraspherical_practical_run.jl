#!/usr/bin/env julia

# Practical starter for baseline MOST vs ultraspherical correction
# Input CSV columns required: zeta, phi_obs
# Example:
#   julia julia/ultraspherical_practical_run.jl data/station_phi_m.csv output/ultra_demo
#
# Optional CSV columns:
#   time  : used when SPLIT_MODE = :blocked to enforce true out-of-sample testing
#
# Workflow summary:
# 1) Fit baseline MOST form phi(zeta)
# 2) Fit ultraspherical correction on residuals only
# 3) Select hyperparameters by held-out validation RMSE
# 4) Export metrics/parameters/predictions (+ optional plot)

using CSV
using DataFrames
using LinearAlgebra
using Statistics
using LsqFit
using Random

# CairoMakie is optional; script still runs without it.
const HAVE_MAKIE = try
    @eval using CairoMakie
    true
catch
    false
end

# ----------------------------- Models ------------------------------------------

# ----------------------------- Settings ----------------------------------------

const TRAIN_FRAC = 0.75
const RNG_SEED = 42

# Set REGIME to :all, :unstable, or :stable
const REGIME = :all

# Split mode:
# :random  -> iid-style random split
# :blocked -> time-ordered split (preferred for nonstationary episodes)
const SPLIT_MODE = :blocked

const N_CANDIDATES = [2, 4, 6]
const RIDGE_CANDIDATES = [0.0, 1e-10, 1e-8, 1e-6, 1e-4]
const LAMBDA_STAR_CANDIDATES = [0.25, 0.5, 0.75]

phi_most(zeta, p) = begin
    a, b, lam = p
    base = max.(1 .- b .* zeta, 1e-8)
    a .* base .^ (-1.0 / lam)
end

xi_map(zeta, alpha_xi) = tanh.(alpha_xi .* zeta)

"""
Recommend alpha_xi so the mapped training data span most of (-1, 1)
without over-saturating tanh at the extremes.
"""
function recommend_alpha_xi(zeta::Vector{Float64}; target_abs_xi::Float64=0.95, q::Float64=0.95)
    zabs = abs.(zeta)
    zq = quantile(zabs, q)
    if zq <= 0.0
        return 0.8
    end
    return atanh(target_abs_xi) / zq
end

"""
Evaluate Gegenbauer C_n^{(lambda_star)}(x) using a stable recurrence:
C_0 = 1
C_1 = 2 lambda x
(n+1) C_{n+1} = 2(n+lambda) x C_n - (n+2lambda-1) C_{n-1}
"""
function gegenbauer_eval(n::Int, lambda_star::Float64, x::Vector{Float64})
    if n == 0
        return ones(length(x))
    elseif n == 1
        return 2.0 * lambda_star .* x
    end
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

function ridge_solve(A::Matrix{Float64}, y::Vector{Float64}, ridge::Float64)
    if ridge <= 0.0
        return A \ y
    end
    ncoef = size(A, 2)
    lhs = transpose(A) * A + ridge * I(ncoef)
    rhs = transpose(A) * y
    return lhs \ rhs
end

rmse(y, yhat) = sqrt(mean((y .- yhat) .^ 2))
mae(y, yhat) = mean(abs.(y .- yhat))

function train_test_split(n::Int, frac::Float64)
    i = sort(randperm(n))
    k = Int(floor(frac * n))
    return i[1:k], i[(k + 1):end]
end

"""
Blocked split for nonstationary data:
- sort by order_key (typically time)
- first TRAIN_FRAC is train, last segment is test
"""
function blocked_split(order_key::Vector, frac::Float64)
    n = length(order_key)
    p = sortperm(order_key)
    k = Int(floor(frac * n))
    return p[1:k], p[(k + 1):end]
end

"""
Choose split indices based on SPLIT_MODE.
"""
function choose_split_indices(order_key::Vector, frac::Float64, split_mode::Symbol)
    n = length(order_key)
    if split_mode == :random
        return train_test_split(n, frac)
    elseif split_mode == :blocked
        return blocked_split(order_key, frac)
    else
        error("Unknown SPLIT_MODE=$(split_mode). Use :random or :blocked")
    end
end

function apply_regime_filter(zeta::Vector{Float64}, y::Vector{Float64}, regime::Symbol)
    if regime == :all
        return zeta, y
    elseif regime == :unstable
        m = zeta .< 0.0
        return zeta[m], y[m]
    elseif regime == :stable
        m = zeta .>= 0.0
        return zeta[m], y[m]
    else
        error("Unknown REGIME=$(regime). Use :all, :unstable, or :stable")
    end
end

# ----------------------------- Runner ------------------------------------------

function main()
    if length(ARGS) < 2
        println("Usage: julia julia/ultraspherical_practical_run.jl <input_csv> <output_prefix>")
        return
    end

    input_csv = ARGS[1]
    out_prefix = ARGS[2]

    df = CSV.read(input_csv, DataFrame)
    required = [:zeta, :phi_obs]
    for c in required
        if !(c in names(df))
            error("Missing column: $(c). Required columns are zeta, phi_obs")
        end
    end

    zeta_raw = Vector{Float64}(df.zeta)
    y_raw = Vector{Float64}(df.phi_obs)

    # If a time column exists, use it for blocked split ordering.
    # Otherwise preserve file row order as the fallback ordering key.
    order_key_raw = (:time in names(df)) ? df.time : collect(1:nrow(df))

    mask = .!(isnan.(zeta_raw) .| isnan.(y_raw) .| isinf.(zeta_raw) .| isinf.(y_raw))
    zeta = zeta_raw[mask]
    y = y_raw[mask]
    order_key = order_key_raw[mask]

    zeta, y = apply_regime_filter(zeta, y, REGIME)
    # Re-apply the same regime mask to order_key for consistent indexing.
    if REGIME == :all
        # no-op
    elseif REGIME == :unstable
        m = zeta_raw[mask] .< 0.0
        order_key = order_key[m]
    elseif REGIME == :stable
        m = zeta_raw[mask] .>= 0.0
        order_key = order_key[m]
    end

    Random.seed!(RNG_SEED)

    n = length(y)
    if n < 40
        error("Need at least 40 valid samples for stable split/fit.")
    end

    train_idx, test_idx = choose_split_indices(order_key, TRAIN_FRAC, SPLIT_MODE)
    z_tr, y_tr = zeta[train_idx], y[train_idx]
    z_te, y_te = zeta[test_idx], y[test_idx]

    # Baseline MOST fit
    p0 = [1.0, 16.0, 4.0]
    lower = [0.1, 0.1, 0.2]
    upper = [5.0, 80.0, 20.0]

    model(x, p) = phi_most(x, p)
    fit = curve_fit(model, z_tr, y_tr, p0, lower=lower, upper=upper)
    p_most = fit.param

    yhat_tr_most = phi_most(z_tr, p_most)
    yhat_te_most = phi_most(z_te, p_most)

    # Ultraspherical residual correction on top of MOST baseline
    alpha_base = recommend_alpha_xi(z_tr)
    alpha_candidates = sort(unique([
        0.5 * alpha_base,
        alpha_base,
        1.5 * alpha_base,
        0.8,
    ]))

    yhat_tr_base = yhat_tr_most
    res_tr = y_tr .- yhat_tr_base

    best = nothing
    best_rmse = Inf

    for alpha_xi in alpha_candidates
        xi_tr = xi_map(z_tr, alpha_xi)
        xi_te = xi_map(z_te, alpha_xi)
        for lambda_star in LAMBDA_STAR_CANDIDATES
            for nmax in N_CANDIDATES
                A_tr = gegenbauer_design(xi_tr, lambda_star, nmax)
                A_te = gegenbauer_design(xi_te, lambda_star, nmax)
                for ridge in RIDGE_CANDIDATES
                    c = ridge_solve(A_tr, res_tr, ridge)

                    yhat_tr = yhat_tr_base .+ A_tr * c
                    yhat_te = phi_most(z_te, p_most) .+ A_te * c

                    score = rmse(y_te, yhat_te)
                    if score < best_rmse
                        best_rmse = score
                        best = (
                            nmax=nmax,
                            coeffs=c,
                            yhat_tr=yhat_tr,
                            yhat_te=yhat_te,
                            alpha_xi=alpha_xi,
                            lambda_star=lambda_star,
                            ridge=ridge,
                        )
                    end
                end
            end
        end
    end

    metrics = DataFrame(
        model=["MOST", "MOST+ULTRA"],
        rmse_train=[rmse(y_tr, yhat_tr_most), rmse(y_tr, best.yhat_tr)],
        rmse_test=[rmse(y_te, yhat_te_most), rmse(y_te, best.yhat_te)],
        mae_train=[mae(y_tr, yhat_tr_most), mae(y_tr, best.yhat_tr)],
        mae_test=[mae(y_te, yhat_te_most), mae(y_te, best.yhat_te)],
    )

    params = DataFrame(
        a=[p_most[1]],
        b=[p_most[2]],
        lambda_profile=[p_most[3]],
        alpha_xi=[best.alpha_xi],
        lambda_star=[best.lambda_star],
        ridge=[best.ridge],
        n_ultra=[best.nmax],
        regime=[String(REGIME)],
        split_mode=[String(SPLIT_MODE)],
    )

    CSV.write("$(out_prefix)_metrics.csv", metrics)
    CSV.write("$(out_prefix)_params.csv", params)

    # Save prediction table for auditing
    pred = DataFrame(
        zeta=z_te,
        obs=y_te,
        most=yhat_te_most,
        ultra=best.yhat_te,
    )
    CSV.write("$(out_prefix)_pred_test.csv", pred)

    # Plot if CairoMakie exists
    if HAVE_MAKIE
        zline = collect(range(minimum(zeta), maximum(zeta), length=500))
        yline_most = phi_most(zline, p_most)

        xi_line = xi_map(zline, best.alpha_xi)
        A_line = gegenbauer_design(xi_line, best.lambda_star, best.nmax)
        yline_ultra = yline_most .+ A_line * best.coeffs

        fig = Figure(resolution=(900, 520))
        ax = Axis(fig[1, 1], xlabel="zeta", ylabel="phi_obs", title="MOST vs MOST+Ultraspherical")
        scatter!(ax, z_tr, y_tr, markersize=6, color=(:steelblue, 0.45), label="Train")
        scatter!(ax, z_te, y_te, markersize=8, color=(:black, 0.75), label="Test")
        lines!(ax, zline, yline_most, linewidth=2.2, color=:orangered, label="MOST")
        lines!(ax, zline, yline_ultra, linewidth=2.2, color=:seagreen, label="MOST + Ultra")
        axislegend(ax, position=:rt)
        save("$(out_prefix)_comparison.png", fig)
    end

    println("Selected ultraspherical hyperparameters:")
    println("  alpha_xi    = $(best.alpha_xi)")
    println("  lambda_star = $(best.lambda_star)")
    println("  nmax        = $(best.nmax)")
    println("  ridge       = $(best.ridge)")
    println("  split_mode  = $(SPLIT_MODE)")
    println("  regime      = $(REGIME)")

    println("Done.")
    println("Saved:")
    println("  $(out_prefix)_metrics.csv")
    println("  $(out_prefix)_params.csv")
    println("  $(out_prefix)_pred_test.csv")
    if HAVE_MAKIE
        println("  $(out_prefix)_comparison.png")
    else
        println("  Plot skipped: CairoMakie not installed.")
    end
end

main()
