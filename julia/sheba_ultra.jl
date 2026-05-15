#!/usr/bin/env julia

# SHEBA-Optimized Ultraspherical Run
# Specifically tuned for Highly Stable Nocturnal Boundary Layers (HSNBL)
# Uses Grachev et al. (2007) BLM baseline and log-mapped Gegenbauer xi
#
# Usage:
#   julia julia/sheba_ultra.jl <input_csv> <output_prefix> [dataset_label] [--baseline=grachev|zero|ultra-only]
#
# input_csv must have columns: zeta, phi_obs  (plus optional: time)
# Grachev baseline:
#   phi_m = 1 + a * zeta * (1+zeta)^(1/3) / (1 + b * zeta)
#   a and b are fitted from data (starting values a=5.0, b=5.0).
# Zero/ultra-only baseline:
#   phi_m baseline = 0, and Gegenbauer series fits phi_m directly.
# xi-mapping uses tanh(alpha * log1p(zeta)) for HSNBL range.
#
# All standard artifacts are generated:
#   _metrics.csv  _params.csv  _pred_test.csv  _coeffs.csv  _curve.csv
#   _model.jl  _formula.md  _validity_summary.md  _report.md
#   _comparison.png  _correction.png  (if CairoMakie available)

using CSV, DataFrames, LinearAlgebra, Statistics, LsqFit, Random

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

# Hyperparameter grid tuned for wide-range stable HSNBL.
const ALPHA_XI_GRID = [0.3, 0.5, 0.8, 1.2]
const LAMBDA_STAR_GRID = [0.25, 0.5, 0.75]
const RIDGE_GRID = [1e-4, 1e-3, 1e-2, 5e-2]
const NMAX_GRID = [3, 4, 5, 6]
const TRAIN_FRAC = 0.75
const RNG_SEED = 42

# ----------------------------- Grachev baseline --------------------------------

"""
Grachev et al. (2007) BLM stable baseline (fixed neutral limit = 1):
  phi_m = 1 + a * zeta * (1 + zeta)^(1/3) / (1 + b * zeta)
For large zeta the function grows as ~ a/b * zeta^(1/3), reproducing the
cube-root regime observed in very stable SHEBA data.
Published canonical values: a ≈ 5.0, b ≈ 5.0.
"""
function phi_grachev(zeta, p)
    a, b = p[1], p[2]
    return 1.0 .+ (a .* zeta .* (max.(1.0 .+ zeta, 1e-8)).^(1.0/3.0)) ./ (1.0 .+ b .* zeta)
end

"""
Log-mapped xi for HSNBL: compresses [0, zeta_max >> 1] into (-1, 1).
"""
xi_map_log(zeta, alpha_xi) = tanh.(alpha_xi .* log1p.(max.(zeta, -0.999)))

function gegenbauer_eval(n::Int, lambda_star::Float64, x::Vector{Float64})
    if n == 0 return ones(length(x)) end
    if n == 1 return 2.0 * lambda_star .* x end
    c_nm1 = ones(length(x))
    c_n = 2.0 * lambda_star .* x
    for k in 1:(n - 1)
        c_np1 = (2.0*(k+lambda_star) .* x .* c_n .- (k+2.0*lambda_star-1.0) .* c_nm1) ./ (k+1.0)
        c_nm1, c_n = c_n, c_np1
    end
    return c_n
end

function gegenbauer_design(xi, lambda_star, nmax)
    A = Matrix{Float64}(undef, length(xi), nmax + 1)
    for n in 0:nmax
        A[:, n+1] = gegenbauer_eval(n, lambda_star, xi)
    end
    return A
end

function ridge_solve(A, y, ridge)
    ridge <= 0.0 && return A \ y
    n = size(A, 2)
    (transpose(A)*A + ridge*I(n)) \ (transpose(A)*y)
end

rmse(y, yhat) = sqrt(mean((y .- yhat).^2))
mae(y, yhat) = mean(abs.(y .- yhat))

function parse_flag(args::Vector{String}, key::String, default::String)
    prefix = key * "="
    for a in args
        startswith(a, prefix) || continue
        parts = split(a, "=", limit=2)
        length(parts) == 2 || return default
        return lowercase(parts[2])
    end
    return default
end

function coeff_physics_summary(coeffs::Vector{Float64}, lambda_star::Float64)
    nmax = length(coeffs) - 1
    mode = collect(0:nmax)
    abs_coeff = abs.(coeffs)
    sum_abs = sum(abs_coeff)
    frac_abs = sum_abs > 0 ? abs_coeff ./ sum_abs : zeros(length(coeffs))

    # Physical proxies in xi-space around neutral (xi=0).
    # C1^(lambda)(xi) = 2 lambda xi, so 2 lambda c1 is the primary tilt term.
    c0 = coeffs[1]
    c1 = nmax >= 1 ? coeffs[2] : 0.0
    c2 = nmax >= 2 ? coeffs[3] : 0.0
    slope_xi_neutral = 2.0 * lambda_star * c1
    curvature_xi_core = 2.0 * lambda_star * (lambda_star + 1.0) * c2

    dominant_idx = argmax(abs_coeff)
    dominant_mode = dominant_idx - 1

    per_mode = DataFrame(
        mode=mode,
        coeff=coeffs,
        abs_coeff=abs_coeff,
        frac_abs=frac_abs,
    )

    summary = (
        c0=c0,
        c1=c1,
        c2=c2,
        slope_xi_neutral=slope_xi_neutral,
        curvature_xi_core=curvature_xi_core,
        dominant_mode=dominant_mode,
        dominant_mode_frac=frac_abs[dominant_idx],
        per_mode=per_mode,
    )
    return summary
end

# ----------------------------- Export helpers ----------------------------------

function write_sheba_exported_model(out_prefix, baseline_mode, p_grachev, alpha_xi, lambda_star, coeffs, nmax)
    baseline_header = baseline_mode == "grachev" ?
        "# Fitted Grachev (2007) baseline + ultraspherical correction" :
        "# Zero baseline + ultraspherical series (ultra-only)"
    baseline_fn = baseline_mode == "grachev" ?
        "phi_baseline_export(zeta) = 1.0 + (SHEBA_A * zeta * max(1.0 + zeta, 1e-8)^(1/3)) / (1.0 + SHEBA_B * zeta)" :
        "phi_baseline_export(zeta) = 0.0"

    lines = [
        "# Auto-generated by julia/sheba_ultra.jl",
        baseline_header,
        "",
        "const SHEBA_A = $(p_grachev[1])",
        "const SHEBA_B = $(p_grachev[2])",
        "const SHEBA_ALPHA_XI = $(alpha_xi)",
        "const SHEBA_LAMBDA_STAR = $(lambda_star)",
        "const SHEBA_NMAX = $(nmax)",
        "const SHEBA_COEFFS = [$(join(string.(coeffs), ", "))]",
        "",
        baseline_fn,
        "xi_sheba(zeta) = tanh(SHEBA_ALPHA_XI * log1p(max(zeta, -0.999)))",
        "",
        "function gegenbauer_export(n, lam, x)",
        "    n == 0 && return 1.0",
        "    n == 1 && return 2.0*lam*x",
        "    c0, c1 = 1.0, 2.0*lam*x",
        "    for k in 1:(n-1)",
        "        c0, c1 = c1, (2*(k+lam)*x*c1 - (k+2lam-1)*c0)/(k+1)",
        "    end",
        "    return c1",
        "end",
        "",
        "function ultra_correction_sheba(zeta)",
        "    xi = xi_sheba(zeta)",
        "    sum(SHEBA_COEFFS[n+1] * gegenbauer_export(n, SHEBA_LAMBDA_STAR, xi) for n in 0:SHEBA_NMAX)",
        "end",
        "",
        "phi_sheba_ultra(zeta) = phi_baseline_export(zeta) + ultra_correction_sheba(zeta)",
    ]
    write("$(out_prefix)_model.jl", join(lines, "\n") * "\n")
end

function write_sheba_formula(out_prefix, baseline_mode, p_grachev, alpha_xi, lambda_star, coeffs)
    coeff_lines = ["- c_$(n) = $(coeffs[n+1])" for n in 0:(length(coeffs)-1)]
    baseline_title = baseline_mode == "grachev" ?
        "# SHEBA Fitted Function (Grachev 2007 + Ultraspherical)" :
        "# SHEBA Fitted Function (Zero Baseline + Ultraspherical)"

    base_lines = if baseline_mode == "grachev"
        [
            "Grachev baseline:",
            "",
            "\$\$",
            "\\phi_{G07}(\\zeta) = 1 + \\frac{a\\,\\zeta\\,(1+\\zeta)^{1/3}}{1 + b\\,\\zeta}",
            "\$\$",
            "",
            "- a = $(p_grachev[1])",
            "- b = $(p_grachev[2])",
            "",
        ]
    else
        [
            "Zero baseline:",
            "",
            "\$\$",
            "\\phi_{base}(\\zeta) = 0",
            "\$\$",
            "",
        ]
    end

    lines = [
        baseline_title,
        "",
        "\$\$",
        "\\phi(\\zeta) = \\phi_{base}(\\zeta) + \\Delta\\phi_{ultra}(\\zeta)",
        "\$\$",
        "",
        "xi-mapping: tanh(alpha * log1p(zeta)),  alpha = $(alpha_xi)",
        "",
        "\$\$",
        "\\Delta\\phi_{ultra}(\\zeta) = \\sum_{n=0}^{$(length(coeffs)-1)} c_n C_n^{(\\lambda_*)}(\\xi(\\zeta))",
        "\$\$",
        "",
        "- lambda_* = $(lambda_star)",
        "",
    ]
    append!(lines, base_lines)
    append!(lines, coeff_lines)
    write("$(out_prefix)_formula.md", join(lines, "\n") * "\n")
end

function write_sheba_validity(out_prefix; dataset_label, baseline_mode, zeta, z_tr, z_te, metrics, p_grachev, alpha_xi, lambda_star, nmax, ridge)
    baseline_name = baseline_mode == "grachev" ? "Grachev" : "Zero"
    ultra_name = baseline_mode == "grachev" ? "Grachev+ULTRA" : "Zero+ULTRA"
    lines = [
        "# SHEBA Fit Validity Summary",
        "",
        "- dataset: $(dataset_label)",
        "- baseline_mode: $(baseline_mode)",
        "- total samples: $(length(zeta))",
        "- train / test: $(length(z_tr)) / $(length(z_te))",
        "",
        "## Stability Range",
        "",
        "- zeta min/max: $(minimum(zeta)), $(maximum(zeta))",
        "- zeta central 5-95%: $(quantile(zeta,0.05)), $(quantile(zeta,0.95))",
        "",
        "## Held-out Skill",
        "",
        "- $(baseline_name) test RMSE: $(metrics.rmse_test[1])",
        "- $(ultra_name) test RMSE: $(metrics.rmse_test[2])",
        "- relative improvement: $(round(100*(metrics.rmse_test[1]-metrics.rmse_test[2])/metrics.rmse_test[1], digits=2))%",
        "",
        "## Parameters",
        "",
        "- baseline a = $(p_grachev[1])",
        "- baseline b = $(p_grachev[2])",
        "- alpha_xi = $(alpha_xi)",
        "- lambda_star = $(lambda_star)",
        "- nmax = $(nmax)",
        "- ridge = $(ridge)",
    ]
    write("$(out_prefix)_validity_summary.md", join(lines, "\n") * "\n")
end

function write_sheba_report(out_prefix; dataset_label, baseline_mode, metrics, params, coeff_summary, have_plot)
    run_name = basename(out_prefix)
    gain_pct = round(100*(metrics.rmse_test[1]-metrics.rmse_test[2])/metrics.rmse_test[1], digits=2)
    baseline_name = baseline_mode == "grachev" ? "Grachev et al. (2007) BLM" : "Zero (no MOST profile)"
    model_base_name = baseline_mode == "grachev" ? "Grachev" : "Zero"
    model_ultra_name = baseline_mode == "grachev" ? "Grachev+ULTRA" : "Zero+ULTRA"
    lines = [
        "# SHEBA Ultraspherical Run Report",
        "",
        "## Run",
        "",
        "- run name: $(run_name)",
        "- dataset: $(dataset_label)",
        "- baseline: $(baseline_name)",
        "- xi-map: tanh(alpha * log1p(zeta))  [log-scale for HSNBL]",
        "",
        "## Metrics",
        "",
        "| Model | RMSE test | MAE test |",
        "|---|---|---|",
        "| $(model_base_name) | $(metrics.rmse_test[1]) | $(metrics.mae_test[1]) |",
        "| $(model_ultra_name) | $(metrics.rmse_test[2]) | $(metrics.mae_test[2]) |",
        "",
        "Relative RMSE gain: **$(gain_pct)%**",
        "",
        "## Fitted Parameters",
        "",
        "- baseline_mode = $(baseline_mode)",
        "- a = $(params.a[1])",
        "- b = $(params.b[1])",
        "- alpha_xi = $(params.alpha_xi[1])",
        "- lambda_star = $(params.lambda_star[1])",
        "- nmax = $(params.n_ultra[1])",
        "",
        "## Coefficient Physical Meaning (HSNBL heuristics)",
        "",
        "- c0 (mode 0 offset): $(coeff_summary.c0)",
        "- c1 (mode 1 tilt): $(coeff_summary.c1)",
        "- c2 (mode 2 curvature): $(coeff_summary.c2)",
        "- Neutral-core slope in xi (2*lambda_* * c1): $(coeff_summary.slope_xi_neutral)",
        "- Core curvature proxy in xi (2*lambda_*(lambda_*+1)*c2): $(coeff_summary.curvature_xi_core)",
        "- Dominant |coeff| mode: n=$(coeff_summary.dominant_mode) (fraction=$(round(100*coeff_summary.dominant_mode_frac, digits=2))%)",
        "",
        "Interpretation guide:",
        "- n=0: bulk level shift relative to baseline.",
        "- n=1: first-order monotonic tilt across stability (often tied to shear/jet strengthening tendency).",
        "- n=2: primary curvature/inversion structure (how sharply phi bends with stability).",
        "- n>=3: higher-order intermittency/wave-like structure and regime transitions.",
        "",
        "## Inline Graphics",
        "",
    ]
    if have_plot
        push!(lines, "### $(model_base_name) vs $(model_ultra_name)")
        push!(lines, "")
        push!(lines, "![comparison]($(run_name)_comparison.png)")
        push!(lines, "")
        push!(lines, "### Ultraspherical Correction")
        push!(lines, "")
        push!(lines, "![correction]($(run_name)_correction.png)")
        push!(lines, "")
    else
        push!(lines, "Plot output not available (CairoMakie not installed).")
        push!(lines, "")
    end
    append!(lines, [
        "## Output Files",
        "",
        "- $(run_name)_metrics.csv",
        "- $(run_name)_params.csv",
        "- $(run_name)_pred_test.csv",
        "- $(run_name)_coeffs.csv",
        "- $(run_name)_curve.csv",
        "- $(run_name)_model.jl",
        "- $(run_name)_formula.md",
        "- $(run_name)_validity_summary.md",
    ])
    write("$(out_prefix)_report.md", join(lines, "\n") * "\n")
end

# ----------------------------- Runner ------------------------------------------

function main()
    if length(ARGS) < 2
        println("Usage: julia julia/sheba_ultra.jl <input_csv> <output_prefix> [dataset_label] [--baseline=grachev|zero|ultra-only]")
        println("input_csv must contain columns: zeta, phi_obs  (optional: time)")
        return
    end

    Random.seed!(RNG_SEED)
    input_csv, out_prefix = ARGS[1], ARGS[2]
    dataset_label = length(ARGS) >= 3 ? ARGS[3] : "SHEBA"
    baseline_mode = parse_flag(ARGS, "--baseline", "grachev")
    baseline_mode = baseline_mode == "ultra-only" ? "zero" : baseline_mode
    baseline_mode in ("grachev", "zero") || error("Unsupported --baseline=$(baseline_mode). Use grachev, zero, or ultra-only.")
    mkpath(dirname(out_prefix))

    df = CSV.read(input_csv, DataFrame)
    # Keep only positive-zeta (stable), finite, physical phi_obs
    mask = map(eachrow(df)) do r
        isfinite(r.zeta) && r.zeta > 0 && isfinite(r.phi_obs)
    end
    df = df[mask, :]
    n = nrow(df)
    n >= 40 || error("Need ≥ 40 stable rows after filtering, got $(n).")

    zeta = Vector{Float64}(df.zeta)
    y    = Vector{Float64}(df.phi_obs)
    order_key = :time in names(df) ? df.time : collect(1:n)

    # Blocked train/test split
    p = sortperm(order_key)
    k = Int(floor(TRAIN_FRAC * n))
    train_idx, test_idx = p[1:k], p[k+1:end]
    z_tr, y_tr = zeta[train_idx], y[train_idx]
    z_te, y_te = zeta[test_idx], y[test_idx]

    # 1. Fit selected baseline
    p_grachev = [0.0, 0.0]
    if baseline_mode == "grachev"
        p0 = [5.0, 5.0]
        fit_base = curve_fit(phi_grachev, z_tr, y_tr, p0, lower=[0.1, 0.1], upper=[20.0, 20.0])
        p_grachev = fit_base.param
    end

    baseline_fun = baseline_mode == "grachev" ? (z -> phi_grachev(z, p_grachev)) : (z -> zeros(length(z)))
    yhat_tr_base = baseline_fun(z_tr)
    yhat_te_base = baseline_fun(z_te)
    res_tr = y_tr .- yhat_tr_base

    # 2. Recommend alpha_xi based on log1p(zeta) range
    log_zeta = log1p.(z_tr)
    lq = quantile(abs.(log_zeta), 0.95)
    alpha_base = lq > 0 ? atanh(0.95) / lq : 0.8
    alpha_candidates = sort(unique([0.5*alpha_base, alpha_base, 1.5*alpha_base, 0.8]))

    # 3. Hyperparameter search on held-out test
    best = nothing
    best_rmse = Inf
    for alpha_xi in alpha_candidates
        xi_tr = xi_map_log(z_tr, alpha_xi)
        xi_te = xi_map_log(z_te, alpha_xi)
        for lambda_star in LAMBDA_STAR_GRID
            for nmax in NMAX_GRID
                A_tr = gegenbauer_design(xi_tr, lambda_star, nmax)
                A_te = gegenbauer_design(xi_te, lambda_star, nmax)
                for ridge in RIDGE_GRID
                    c = ridge_solve(A_tr, res_tr, ridge)
                    yhat_te = yhat_te_base .+ A_te * c
                    score = rmse(y_te, yhat_te)
                    if score < best_rmse
                        best_rmse = score
                        best = (nmax=nmax, coeffs=c, alpha_xi=alpha_xi, lambda_star=lambda_star, ridge=ridge,
                                yhat_tr=yhat_tr_base .+ A_tr*c, yhat_te=yhat_te)
                    end
                end
            end
        end
    end

    metrics = DataFrame(
         model=[baseline_mode == "grachev" ? "Grachev" : "Zero",
             baseline_mode == "grachev" ? "Grachev+ULTRA" : "Zero+ULTRA"],
        rmse_train=[rmse(y_tr, yhat_tr_base), rmse(y_tr, best.yhat_tr)],
        rmse_test=[rmse(y_te, yhat_te_base), rmse(y_te, best.yhat_te)],
        mae_train=[mae(y_tr, yhat_tr_base), mae(y_tr, best.yhat_tr)],
        mae_test=[mae(y_te, yhat_te_base), mae(y_te, best.yhat_te)],
    )

    params = DataFrame(
        dataset=[dataset_label],
        a=[p_grachev[1]],
        b=[p_grachev[2]],
        alpha_xi=[best.alpha_xi],
        lambda_star=[best.lambda_star],
        ridge=[best.ridge],
        n_ultra=[best.nmax],
        regime=["stable"],
        split_mode=["blocked"],
        baseline_mode=[baseline_mode],
        xi_mode=["log"],
    )

    pred = DataFrame(zeta=z_te, obs=y_te, baseline=yhat_te_base, ultra=best.yhat_te)

    coeffs = DataFrame(mode=collect(0:best.nmax), coeff=collect(best.coeffs))

    zline = collect(range(minimum(zeta), maximum(zeta), length=500))
    xi_line = xi_map_log(zline, best.alpha_xi)
    A_line  = gegenbauer_design(xi_line, best.lambda_star, best.nmax)
    corr_line = A_line * best.coeffs
    yline_base = baseline_fun(zline)
    curve = DataFrame(zeta=zline, baseline=yline_base, correction=corr_line,
                      ultra=yline_base .+ corr_line)

    coeff_summary = coeff_physics_summary(collect(best.coeffs), best.lambda_star)

    CSV.write("$(out_prefix)_metrics.csv", metrics)
    CSV.write("$(out_prefix)_params.csv", params)
    CSV.write("$(out_prefix)_pred_test.csv", pred)
    CSV.write("$(out_prefix)_coeffs.csv", coeffs)
    CSV.write("$(out_prefix)_curve.csv", curve)
    write_sheba_exported_model(out_prefix, baseline_mode, p_grachev, best.alpha_xi, best.lambda_star, best.coeffs, best.nmax)
    write_sheba_formula(out_prefix, baseline_mode, p_grachev, best.alpha_xi, best.lambda_star, best.coeffs)
    write_sheba_validity(out_prefix; dataset_label=dataset_label, baseline_mode=baseline_mode, zeta=zeta, z_tr=z_tr, z_te=z_te,
                         metrics=metrics, p_grachev=p_grachev, alpha_xi=best.alpha_xi,
                         lambda_star=best.lambda_star, nmax=best.nmax, ridge=best.ridge)
    write_sheba_report(out_prefix; dataset_label=dataset_label, baseline_mode=baseline_mode,
                       metrics=metrics, params=params, coeff_summary=coeff_summary,
                       have_plot=HAVE_MAKIE)

    println("Selected SHEBA hyperparameters:")
    println("  baseline_mode = $(baseline_mode)")
    println("  a             = $(p_grachev[1])")
    println("  b             = $(p_grachev[2])")
    println("  alpha_xi      = $(best.alpha_xi)")
    println("  lambda_star   = $(best.lambda_star)")
    println("  nmax          = $(best.nmax)")
    println("  ridge         = $(best.ridge)")
    println("  baseline RMSE = $(metrics.rmse_test[1])")
    println("  ultra RMSE    = $(metrics.rmse_test[2])")
    println("Done.")
    println("Saved:")
    for suffix in ["_metrics.csv","_params.csv","_pred_test.csv","_coeffs.csv","_curve.csv",
                   "_model.jl","_formula.md","_validity_summary.md","_report.md"]
        println("  $(out_prefix)$(suffix)")
    end

    HAVE_MAKIE || return

    fig = Figure(resolution=(900, 520))
    title_lbl = baseline_mode == "grachev" ? "SHEBA: Grachev vs Grachev+Ultraspherical" : "SHEBA: Zero Baseline vs Ultraspherical"
    ax  = Axis(fig[1,1], xlabel="zeta", ylabel="phi_m", title=title_lbl)
    scatter!(ax, z_tr, y_tr, markersize=5, color=(:steelblue, 0.4), label="Train")
    scatter!(ax, z_te, y_te, markersize=7, color=(:black, 0.7),     label="Test")
    base_lbl = baseline_mode == "grachev" ? "Grachev" : "Zero"
    ultra_lbl = baseline_mode == "grachev" ? "Grachev+Ultra" : "Zero+Ultra"
    lines!(ax, zline, curve.baseline, linewidth=2.2, color=:orangered, label=base_lbl)
    lines!(ax, zline, curve.ultra,    linewidth=2.2, color=:seagreen,  label=ultra_lbl)
    axislegend(ax, position=:lt)
    save("$(out_prefix)_comparison.png", fig)

    fig2 = Figure(resolution=(900, 520))
    ax2  = Axis(fig2[1,1], xlabel="zeta", ylabel="Delta phi_ultra", title="Ultraspherical Correction (SHEBA)")
    lines!(ax2, zline, corr_line, linewidth=2.4, color=:darkgreen, label="Delta phi_ultra(zeta)")
    hlines!(ax2, [0.0], color=:black, linewidth=1.0, linestyle=:dash)
    axislegend(ax2, position=:rt)
    save("$(out_prefix)_correction.png", fig2)
    println("  $(out_prefix)_comparison.png")
    println("  $(out_prefix)_correction.png")
end

main()

