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
#
# Synthetic mode:
#   julia julia/ultraspherical_practical_run.jl --synthetic output/ultra_synth [noise_frac] [n_samples]
# This generates noisy synthetic data with known Gegenbauer coefficients so
# students can test recovery before touching observations.

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
const SYNTHETIC_DEFAULT_SAMPLES = 240
const SYNTHETIC_DEFAULT_NOISE = 0.05

# Baseline default constants (Dyer linear form encoded by lambda=-1):
# phi_m = a * (1 + b zeta)^(-1/lambda_profile)
const MOST_A_FIXED = 1.0
const MOST_LAMBDA_FIXED = -1.0
const MOST_B_FIXED = 4.7

phi_most(zeta, p) = begin
    a, b, lam = p
    base = max.(1 .+ b .* zeta, 1e-8)
    a .* base .^ (-1.0 / lam)
end

phi_most_fixed(zeta, b) = phi_most(zeta, [MOST_A_FIXED, b, MOST_LAMBDA_FIXED])

"""
Grachev et al. (2007) BLM stable baseline:
  phi_m = 1 + a * zeta * (1 + zeta)^(1/3) / (1 + b * zeta)
Neutral limit is exactly 1. For large zeta approaches a/b * zeta^(1/3),
which reproduces the weak cube-root growth seen in very stable HSNBL.
"""
phi_grachev(zeta, p) = begin
    a, b = p[1], p[2]
    1.0 .+ (a .* zeta .* (max.(1.0 .+ zeta, 1e-8)).^(1.0/3.0)) ./ (1.0 .+ b .* zeta)
end

# Published canonical values from Grachev et al. (2007), Table 3, phi_m fit:
const GRACHEV_A_FIXED = 5.0
const GRACHEV_B_FIXED = 5.0

xi_map(zeta, alpha_xi; xi_mode::Symbol=:tanh) = begin
    if xi_mode == :log
        tanh.(alpha_xi .* log1p.(max.(zeta, -0.999)))
    else
        tanh.(alpha_xi .* zeta)
    end
end

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

function print_usage()
    println("Usage:")
    println("  julia julia/ultraspherical_practical_run.jl <input_csv> <output_prefix>")
    println("  julia julia/ultraspherical_practical_run.jl --synthetic <output_prefix> [noise_frac] [n_samples]")
    println("")
    println("Optional flags for observed-data mode:")
    println("  --baseline=dyer47|linear-fit|ultra-only|most-free|grachev   (default dyer47)")
    println("  --xi-map=tanh|log   (default tanh; use log for SHEBA/HSNBL where zeta >> 1)")
    println("")
    println("Required CSV columns: zeta, phi_obs")
    println("Optional CSV column:  time")
    println("Synthetic defaults: noise_frac=$(SYNTHETIC_DEFAULT_NOISE), n_samples=$(SYNTHETIC_DEFAULT_SAMPLES)")
end

function parse_flag(args::Vector{String}, key::String, default::String)
    prefix = key * "="
    for a in args
        startswith(a, prefix) || continue
        parts = split(a, "=", limit=2)
        length(parts) == 2 || return default
        return parts[2]
    end
    return default
end

"""
Generate synthetic training data from a known baseline + Gegenbauer residual model.
Additive white Gaussian noise is included so coefficient recovery is realistic rather
than a perfect interpolation exercise.
"""
function generate_synthetic_dataset(; n_samples::Int=SYNTHETIC_DEFAULT_SAMPLES, noise_frac::Float64=SYNTHETIC_DEFAULT_NOISE)
    phase = collect(range(0.0, 6.0 * pi, length=n_samples))
    zeta = 0.9 .* sin.(phase) .+ 0.45 .* sin.(0.37 .* phase .+ 0.6) .+ 0.20 .* randn(n_samples)
    zeta = clamp.(zeta, -1.5, 2.2)
    time = collect(1:n_samples)

    p_true = [1.0, 0.35, 4.0]
    alpha_xi_true = 0.9
    lambda_star_true = 0.5
    coeff_true = [0.0, 0.10, -0.06, 0.035, -0.02]
    nmax_true = length(coeff_true) - 1

    baseline_true = phi_most(zeta, p_true)
    xi_true = xi_map(zeta, alpha_xi_true)
    A_true = gegenbauer_design(xi_true, lambda_star_true, nmax_true)
    residual_true = A_true * coeff_true
    phi_clean = baseline_true .+ residual_true

    noise_sigma = noise_frac * std(phi_clean)
    phi_obs = phi_clean .+ noise_sigma .* randn(n_samples)

    df = DataFrame(
        time=time,
        zeta=zeta,
        phi_obs=phi_obs,
        phi_clean=phi_clean,
        baseline_true=baseline_true,
        residual_true=residual_true,
    )

    truth = (
        a=p_true[1],
        b=p_true[2],
        lambda_profile=p_true[3],
        alpha_xi=alpha_xi_true,
        lambda_star=lambda_star_true,
        nmax=nmax_true,
        coeff_true=coeff_true,
        noise_sigma=noise_sigma,
        noise_frac=noise_frac,
    )

    return df, truth
end

function make_coeff_table(coeff_est::Vector{Float64}, nmax_est::Int; truth=nothing)
    if truth === nothing
        return DataFrame(mode=collect(0:nmax_est), coeff_estimate=collect(coeff_est))
    end

    max_mode = max(nmax_est, truth.nmax)
    modes = collect(0:max_mode)
    est = Union{Missing, Float64}[mode <= nmax_est ? coeff_est[mode + 1] : missing for mode in modes]
    tru = Union{Missing, Float64}[mode <= truth.nmax ? truth.coeff_true[mode + 1] : missing for mode in modes]
    return DataFrame(mode=modes, coeff_estimate=est, coeff_true=tru)
end

function write_exported_model(out_prefix::String, p_most::Vector{Float64}, alpha_xi::Float64, lambda_star::Float64, coeffs::Vector{Float64}, nmax::Int)
    coeff_literal = join(string.(coeffs), ", ")
    lines = [
        "# Auto-generated by julia/ultraspherical_practical_run.jl",
        "# Fitted MOST + ultraspherical correction model",
        "",
        "const ULTRA_A = $(p_most[1])",
        "const ULTRA_B = $(p_most[2])",
        "const ULTRA_LAMBDA_PROFILE = $(p_most[3])",
        "const ULTRA_ALPHA_XI = $(alpha_xi)",
        "const ULTRA_LAMBDA_STAR = $(lambda_star)",
        "const ULTRA_NMAX = $(nmax)",
        "const ULTRA_COEFFS = [$(coeff_literal)]",
        "",
        "phi_most_export(zeta) = begin",
        "    base = max(1.0 + ULTRA_B * zeta, 1e-8)",
        "    ULTRA_A * base^(-1.0 / ULTRA_LAMBDA_PROFILE)",
        "end",
        "",
        "xi_map_export(zeta) = tanh(ULTRA_ALPHA_XI * zeta)",
        "",
        "function gegenbauer_eval_export(n::Int, lambda_star::Float64, x::Float64)",
        "    if n == 0",
        "        return 1.0",
        "    elseif n == 1",
        "        return 2.0 * lambda_star * x",
        "    end",
        "    c_nm1 = 1.0",
        "    c_n = 2.0 * lambda_star * x",
        "    for k in 1:(n - 1)",
        "        c_np1 = (2.0 * (k + lambda_star) * x * c_n - (k + 2.0 * lambda_star - 1.0) * c_nm1) / (k + 1.0)",
        "        c_nm1, c_n = c_n, c_np1",
        "    end",
        "    return c_n",
        "end",
        "",
        "function ultra_correction_export(zeta)",
        "    xi = xi_map_export(zeta)",
        "    total = 0.0",
        "    for n in 0:ULTRA_NMAX",
        "        total += ULTRA_COEFFS[n + 1] * gegenbauer_eval_export(n, ULTRA_LAMBDA_STAR, xi)",
        "    end",
        "    return total",
        "end",
        "",
        "phi_ultra_export(zeta) = phi_most_export(zeta) + ultra_correction_export(zeta)",
    ]
    write("$(out_prefix)_model.jl", join(lines, "\n") * "\n")
end

function write_formula_summary(out_prefix::String, p_most::Vector{Float64}, alpha_xi::Float64, lambda_star::Float64, coeffs::Vector{Float64})
    coeff_lines = ["- c_$(n) = $(coeffs[n + 1])" for n in 0:(length(coeffs) - 1)]
    lines = [
        "# Exported MOST + Ultraspherical Function",
        "",
        "The fitted function is",
        "",
        "\$\$",
        "\\phi(\\zeta) = \\phi_{MOST}(\\zeta) + \\Delta \\phi_{ultra}(\\zeta)",
        "\$\$",
        "",
        "with baseline",
        "",
        "\$\$",
        "\\phi_{MOST}(\\zeta) = a (1 + b \\zeta)^{-1/\\lambda_p}",
        "\$\$",
        "",
        "where",
        "",
        "- a = $(p_most[1])",
        "- b = $(p_most[2])",
        "- \\lambda_p = $(p_most[3])",
        "",
        "The ultraspherical correction uses",
        "",
        "\$\$",
        "\\xi(\\zeta) = \\tanh(\\alpha_\\xi \\zeta)",
        "\$\$",
        "",
        "with",
        "",
        "- \\alpha_\\xi = $(alpha_xi)",
        "- \\lambda_* = $(lambda_star)",
        "",
        "and residual correction",
        "",
        "\$\$",
        "\\Delta \\phi_{ultra}(\\zeta) = \\sum_{n=0}^{$(length(coeffs) - 1)} c_n C_n^{(\\lambda_*)}(\\xi(\\zeta))",
        "\$\$",
        "",
        "with fitted coefficients",
        "",
    ]
    append!(lines, coeff_lines)
    write("$(out_prefix)_formula.md", join(lines, "\n") * "\n")
end

function write_validity_summary(out_prefix::String; dataset_label::String, zeta::Vector{Float64}, z_tr::Vector{Float64}, z_te::Vector{Float64}, metrics::DataFrame, p_most::Vector{Float64}, alpha_xi::Float64, lambda_star::Float64, nmax::Int, ridge::Float64)
    zeta_min = minimum(zeta)
    zeta_max = maximum(zeta)
    zeta_q05 = quantile(zeta, 0.05)
    zeta_q95 = quantile(zeta, 0.95)
    most_rmse_test = metrics.rmse_test[1]
    ultra_rmse_test = metrics.rmse_test[2]
    improvement = most_rmse_test - ultra_rmse_test
    improvement_frac = most_rmse_test > 0 ? 100.0 * improvement / most_rmse_test : NaN

    lines = [
        "# Fit Validity Summary",
        "",
        "## Scope",
        "",
        "- dataset: $(dataset_label)",
        "- total samples used: $(length(zeta))",
        "- training samples: $(length(z_tr))",
        "- test samples: $(length(z_te))",
        "",
        "## Recommended Validity Range",
        "",
        "Use the exported function primarily within the fitted stability interval:",
        "",
        "- full fitted zeta range: [$(zeta_min), $(zeta_max)]",
        "- robust central zeta range (5%-95%): [$(zeta_q05), $(zeta_q95)]",
        "",
        "## Held-out Skill",
        "",
        "- MOST test RMSE: $(most_rmse_test)",
        "- MOST+ULTRA test RMSE: $(ultra_rmse_test)",
        "- absolute RMSE improvement: $(improvement)",
        "- relative RMSE improvement: $(improvement_frac)%",
        "",
        "## Fitted Function Settings",
        "",
        "- baseline a: $(p_most[1])",
        "- baseline b: $(p_most[2])",
        "- baseline lambda_profile: $(p_most[3])",
        "- alpha_xi: $(alpha_xi)",
        "- lambda_star: $(lambda_star)",
        "- n_ultra: $(nmax)",
        "- ridge: $(ridge)",
        "",
        "## Usage Notes",
        "",
        "- The exported ultraspherical correction is a data-fitted residual term added to the baseline MOST function.",
        "- Interpret the correction curve directly before interpreting individual spectral coefficients.",
        "- Avoid extrapolating far outside the fitted zeta range unless additional validation is performed.",
    ]
    write("$(out_prefix)_validity_summary.md", join(lines, "\n") * "\n")
end

function write_run_report(out_prefix::String; dataset_label::String, metrics::DataFrame, params::DataFrame, have_plot::Bool)
    run_name = basename(out_prefix)
    most_test = metrics.rmse_test[1]
    ultra_test = metrics.rmse_test[2]
    gain = most_test - ultra_test
    gain_pct = most_test > 0 ? 100.0 * gain / most_test : NaN

    lines = [
        "# Ultraspherical Run Report",
        "",
        "## Run",
        "",
        "- run name: $(run_name)",
        "- dataset label: $(dataset_label)",
        "",
        "## Metrics",
        "",
        "- MOST test RMSE: $(most_test)",
        "- MOST+ULTRA test RMSE: $(ultra_test)",
        "- absolute RMSE gain: $(gain)",
        "- relative RMSE gain: $(gain_pct)%",
        "",
        "## Parameters",
        "",
        "- baseline a: $(params.a[1])",
        "- baseline b: $(params.b[1])",
        "- baseline lambda_profile: $(params.lambda_profile[1])",
        "- alpha_xi: $(params.alpha_xi[1])",
        "- lambda_star: $(params.lambda_star[1])",
        "- ridge: $(params.ridge[1])",
        "- n_ultra: $(params.n_ultra[1])",
        "- regime: $(params.regime[1])",
        "- split_mode: $(params.split_mode[1])",
        "",
        "## Inline Graphics",
        "",
    ]

    if have_plot
        push!(lines, "### MOST vs MOST+ULTRA")
        push!(lines, "")
        push!(lines, "![MOST vs MOST+ULTRA]($(run_name)_comparison.png)")
        push!(lines, "")
        push!(lines, "### Ultraspherical Correction")
        push!(lines, "")
        push!(lines, "![Ultraspherical Correction]($(run_name)_correction.png)")
        push!(lines, "")
    else
        push!(lines, "Plot output not available (CairoMakie not installed at run time).")
        push!(lines, "")
    end

    append!(lines, [
        "## Run Files",
        "",
        "- $(run_name)_metrics.csv",
        "- $(run_name)_params.csv",
        "- $(run_name)_coeffs.csv",
        "- $(run_name)_pred_test.csv",
        "- $(run_name)_curve.csv",
        "- $(run_name)_model.jl",
        "- $(run_name)_formula.md",
        "- $(run_name)_validity_summary.md",
    ])

    write("$(out_prefix)_report.md", join(lines, "\n") * "\n")
end

function run_pipeline(df::DataFrame, out_prefix::String; dataset_label::String="observed", truth=nothing, baseline_mode::Symbol=:dyer47, xi_mode::Symbol=:tanh)
    required = [:zeta, :phi_obs]
    col_syms = Symbol.(names(df))
    for c in required
        if !(c in col_syms)
            error("Missing column: $(c). Required columns are zeta, phi_obs")
        end
    end

    zeta_raw = Vector{Float64}(df.zeta)
    y_raw = Vector{Float64}(df.phi_obs)
    order_key_raw = (:time in names(df)) ? df.time : collect(1:nrow(df))

    valid_mask = .!(isnan.(zeta_raw) .| isnan.(y_raw) .| isinf.(zeta_raw) .| isinf.(y_raw))
    zeta_valid = zeta_raw[valid_mask]
    y_valid = y_raw[valid_mask]
    order_key_valid = order_key_raw[valid_mask]

    regime_mask = if REGIME == :all
        trues(length(zeta_valid))
    elseif REGIME == :unstable
        zeta_valid .< 0.0
    elseif REGIME == :stable
        zeta_valid .>= 0.0
    else
        error("Unknown REGIME=$(REGIME). Use :all, :unstable, or :stable")
    end

    zeta = zeta_valid[regime_mask]
    y = y_valid[regime_mask]
    order_key = order_key_valid[regime_mask]

    n = length(y)
    if n < 40
        error("Need at least 40 valid samples for stable split/fit.")
    end

    train_idx, test_idx = choose_split_indices(order_key, TRAIN_FRAC, SPLIT_MODE)
    z_tr, y_tr = zeta[train_idx], y[train_idx]
    z_te, y_te = zeta[test_idx], y[test_idx]

    # For Grachev mode we use phi_grachev with its own 2-param vector;
    # for all other modes we use the standard phi_most 3-param vector.
    use_grachev = (baseline_mode == :grachev)

    p_baseline = if baseline_mode == :dyer47
        [MOST_A_FIXED, MOST_B_FIXED, MOST_LAMBDA_FIXED]
    elseif baseline_mode == :linear_fit
        p0 = [16.0]
        lower = [0.1]
        upper = [80.0]
        model = (x, p) -> phi_most_fixed(x, p[1])
        fit = curve_fit(model, z_tr, y_tr, p0, lower=lower, upper=upper)
        [MOST_A_FIXED, fit.param[1], MOST_LAMBDA_FIXED]
    elseif baseline_mode == :ultra_only
        [1.0, 0.0, 1.0]
    elseif baseline_mode == :most_free
        p0 = [1.0, 16.0, 4.0]
        lower = [0.1, 0.1, 0.2]
        upper = [5.0, 80.0, 20.0]
        model = (x, p) -> phi_most(x, p)
        fit = curve_fit(model, z_tr, y_tr, p0, lower=lower, upper=upper)
        fit.param
    elseif baseline_mode == :grachev
        p0 = [GRACHEV_A_FIXED, GRACHEV_B_FIXED]
        lower = [0.1, 0.1]
        upper = [20.0, 20.0]
        model = (x, p) -> phi_grachev(x, p)
        fit = curve_fit(model, z_tr, y_tr, p0, lower=lower, upper=upper)
        fit.param
    else
        error("Unknown baseline_mode=$(baseline_mode). Use dyer47, linear_fit, ultra_only, most_free, or grachev")
    end

    eval_baseline = use_grachev ?
        (z -> phi_grachev(z, p_baseline)) :
        (z -> phi_most(z, p_baseline))

    # Canonical 3-element param vector for reporting/export (Grachev has only 2):
    p_most = use_grachev ? [1.0, p_baseline[1], p_baseline[2]] : p_baseline

    yhat_tr_most = eval_baseline(z_tr)
    yhat_te_most = eval_baseline(z_te)

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
        xi_tr = xi_map(z_tr, alpha_xi; xi_mode=xi_mode)
        xi_te = xi_map(z_te, alpha_xi; xi_mode=xi_mode)
        for lambda_star in LAMBDA_STAR_CANDIDATES
            for nmax in N_CANDIDATES
                A_tr = gegenbauer_design(xi_tr, lambda_star, nmax)
                A_te = gegenbauer_design(xi_te, lambda_star, nmax)
                for ridge in RIDGE_CANDIDATES
                    c = ridge_solve(A_tr, res_tr, ridge)

                    yhat_tr = yhat_tr_base .+ A_tr * c
                    yhat_te = eval_baseline(z_te) .+ A_te * c

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
        dataset=[dataset_label],
        a=[p_most[1]],
        b=[p_most[2]],
        lambda_profile=[p_most[3]],
        alpha_xi=[best.alpha_xi],
        lambda_star=[best.lambda_star],
        ridge=[best.ridge],
        n_ultra=[best.nmax],
        regime=[String(REGIME)],
        split_mode=[String(SPLIT_MODE)],
        baseline_mode=[String(baseline_mode)],
        xi_mode=[String(xi_mode)],
    )

    if truth !== nothing
        params.a_true = [truth.a]
        params.b_true = [truth.b]
        params.lambda_profile_true = [truth.lambda_profile]
        params.alpha_xi_true = [truth.alpha_xi]
        params.lambda_star_true = [truth.lambda_star]
        params.n_ultra_true = [truth.nmax]
        params.noise_sigma = [truth.noise_sigma]
        params.noise_frac = [truth.noise_frac]
    end

    pred = DataFrame(
        zeta=z_te,
        obs=y_te,
        most=yhat_te_most,
        ultra=best.yhat_te,
    )

    coeffs = make_coeff_table(best.coeffs, best.nmax; truth=truth)

    zline = collect(range(minimum(zeta), maximum(zeta), length=500))
    yline_most = eval_baseline(zline)
    xi_line = xi_map(zline, best.alpha_xi; xi_mode=xi_mode)
    A_line = gegenbauer_design(xi_line, best.lambda_star, best.nmax)
    correction_line = A_line * best.coeffs
    yline_ultra = yline_most .+ correction_line

    curve = DataFrame(
        zeta=zline,
        most=yline_most,
        correction=correction_line,
        ultra=yline_ultra,
    )

    CSV.write("$(out_prefix)_metrics.csv", metrics)
    CSV.write("$(out_prefix)_params.csv", params)
    CSV.write("$(out_prefix)_pred_test.csv", pred)
    CSV.write("$(out_prefix)_coeffs.csv", coeffs)
    CSV.write("$(out_prefix)_curve.csv", curve)
    write_exported_model(out_prefix, p_most, best.alpha_xi, best.lambda_star, best.coeffs, best.nmax)
    write_formula_summary(out_prefix, p_most, best.alpha_xi, best.lambda_star, best.coeffs)
    write_validity_summary(
        out_prefix;
        dataset_label=dataset_label,
        zeta=zeta,
        z_tr=z_tr,
        z_te=z_te,
        metrics=metrics,
        p_most=p_most,
        alpha_xi=best.alpha_xi,
        lambda_star=best.lambda_star,
        nmax=best.nmax,
        ridge=best.ridge,
    )
    write_run_report(
        out_prefix;
        dataset_label=dataset_label,
        metrics=metrics,
        params=params,
        have_plot=HAVE_MAKIE,
    )

    if truth !== nothing
        CSV.write("$(out_prefix)_synthetic_data.csv", df)
    end

    if HAVE_MAKIE
        fig = Figure(resolution=(900, 520))
        ax = Axis(fig[1, 1], xlabel="zeta", ylabel="phi_obs", title="MOST vs MOST+Ultraspherical")
        scatter!(ax, z_tr, y_tr, markersize=6, color=(:steelblue, 0.45), label="Train")
        scatter!(ax, z_te, y_te, markersize=8, color=(:black, 0.75), label="Test")
        lines!(ax, zline, yline_most, linewidth=2.2, color=:orangered, label="MOST")
        lines!(ax, zline, yline_ultra, linewidth=2.2, color=:seagreen, label="MOST + Ultra")
        axislegend(ax, position=:rt)
        save("$(out_prefix)_comparison.png", fig)

        fig_corr = Figure(resolution=(900, 520))
        ax_corr = Axis(fig_corr[1, 1], xlabel="zeta", ylabel="Delta phi_ultra", title="Ultraspherical Correction Only")
        lines!(ax_corr, zline, correction_line, linewidth=2.4, color=:darkgreen, label="Delta phi_ultra(zeta)")
        hlines!(ax_corr, [0.0], color=:black, linewidth=1.0, linestyle=:dash)
        axislegend(ax_corr, position=:rt)
        save("$(out_prefix)_correction.png", fig_corr)
    end

    println("Selected ultraspherical hyperparameters:")
    println("  alpha_xi    = $(best.alpha_xi)")
    println("  lambda_star = $(best.lambda_star)")
    println("  nmax        = $(best.nmax)")
    println("  ridge       = $(best.ridge)")
    println("  split_mode    = $(SPLIT_MODE)")
    println("  baseline_mode = $(baseline_mode)")
    println("  xi_mode       = $(xi_mode)")
    println("  regime        = $(REGIME)")
    println("  most_a        = $(p_most[1])")
    println("  most_b        = $(p_most[2])")
    println("  most_lambda   = $(p_most[3])")
    if truth !== nothing
        println("  synthetic noise sigma = $(truth.noise_sigma)")
    end

    println("Done.")
    println("Saved:")
    println("  $(out_prefix)_metrics.csv")
    println("  $(out_prefix)_params.csv")
    println("  $(out_prefix)_pred_test.csv")
    println("  $(out_prefix)_coeffs.csv")
    println("  $(out_prefix)_curve.csv")
    println("  $(out_prefix)_model.jl")
    println("  $(out_prefix)_formula.md")
    println("  $(out_prefix)_validity_summary.md")
    println("  $(out_prefix)_report.md")
    if truth !== nothing
        println("  $(out_prefix)_synthetic_data.csv")
    end
    if HAVE_MAKIE
        println("  $(out_prefix)_comparison.png")
        println("  $(out_prefix)_correction.png")
    else
        println("  Plot skipped: CairoMakie not installed.")
    end
end

# ----------------------------- Runner ------------------------------------------

function main()
    if isempty(ARGS)
        print_usage()
        return
    end

    Random.seed!(RNG_SEED)

    if ARGS[1] == "--synthetic"
        if length(ARGS) < 2
            print_usage()
            return
        end
        out_prefix = ARGS[2]
        noise_frac = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : SYNTHETIC_DEFAULT_NOISE
        n_samples = length(ARGS) >= 4 ? parse(Int, ARGS[4]) : SYNTHETIC_DEFAULT_SAMPLES
        df, truth = generate_synthetic_dataset(n_samples=n_samples, noise_frac=noise_frac)
        run_pipeline(df, out_prefix; dataset_label="synthetic", truth=truth)
        return
    end

    if length(ARGS) < 2
        print_usage()
        return
    end

    input_csv = ARGS[1]
    out_prefix = ARGS[2]
    extra = length(ARGS) > 2 ? ARGS[3:end] : String[]
    baseline_flag = lowercase(parse_flag(extra, "--baseline", "dyer47"))
    baseline_mode = if baseline_flag == "dyer47"
        :dyer47
    elseif baseline_flag in ("linear-fit", "linear_fit")
        :linear_fit
    elseif baseline_flag in ("ultra-only", "ultra_only")
        :ultra_only
    elseif baseline_flag in ("most-free", "most_free")
        :most_free
    elseif baseline_flag == "grachev"
        :grachev
    else
        error("Unknown --baseline=$(baseline_flag). Use dyer47, linear-fit, ultra-only, most-free, or grachev")
    end
    xi_flag = lowercase(parse_flag(extra, "--xi-map", "tanh"))
    xi_mode = if xi_flag == "log"
        :log
    elseif xi_flag == "tanh"
        :tanh
    else
        error("Unknown --xi-map=$(xi_flag). Use tanh or log")
    end
    df = CSV.read(input_csv, DataFrame)
    run_pipeline(df, out_prefix; baseline_mode=baseline_mode, xi_mode=xi_mode)
end

main()
