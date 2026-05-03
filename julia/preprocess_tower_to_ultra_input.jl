#!/usr/bin/env julia

using CSV, DataFrames, Dates, Downloads, Statistics

const KAPPA::Float64 = 0.4
const G::Float64 = 9.81
const DEFAULT_WTV_EPS::Float64 = 1e-12

function print_usage()
    println("Usage:")
    println("  julia julia/preprocess_tower_to_ultra_input.jl <input_or_station> <output_csv> <z_m> <d_m> [--stable-only] [--phi=phi_m|phi_h] [--mode=raw|two-level|api-smear] [--z1=<m> --z2=<m>]")
    println("")
    println("Required raw flux inputs (with alias support):")
    println("  uw covariance: uw, u_w_cov, u_prime_w_prime")
    println("  vw covariance: vw, v_w_cov, v_prime_w_prime")
    println("  buoyancy flux: wthetav, w_theta_v, wthetav_cov")
    println("  mean virtual temperature: thetav, theta_v, theta_v_mean")
    println("")
    println("Mode raw (default): optional direct gradients for phi terms")
    println("  dU_dz aliases: dudz, dU_dz, dUdz")
    println("  dtheta_dz aliases: dthetadz, dtheta_dz, dthetav_dz")
    println("")
    println("Mode two-level: computes gradients from profile levels")
    println("  Required flags: --mode=two-level --z1=<m> --z2=<m>")
    println("  U at z1 aliases: u1, u_low, U1, wind_low")
    println("  U at z2 aliases: u2, u_high, U2, wind_high")
    println("  Theta at z1 aliases: theta1, theta_low, T1, temp_low")
    println("  Theta at z2 aliases: theta2, theta_high, T2, temp_high")
    println("")
    println("Mode api-smear: fetches data directly from SmartSMEAR API")
    println("  Positional arg #1 is station label for logs (e.g., HYY)")
    println("  Required flags: --mode=api-smear --from=ISO --to=ISO")
    println("  Optional flags: --interval=30 --aggregation=ARITHMETIC --quality=ANY")
    println("  Neutral handling: --wtv-eps=<eps> (default $(DEFAULT_WTV_EPS))")
    println("  Optional direct stability inputs:")
    println("    --tv-mo-length=TABLE.VAR [--tv-ustar=TABLE.VAR]")
    println("    When --tv-mo-length is provided, zeta uses z_eff / MO_length directly.")
    println("  Required tablevariable flags for flux terms:")
    println("    --tv-uw=TABLE.VAR --tv-vw=TABLE.VAR --tv-wthetav=TABLE.VAR --tv-thetav=TABLE.VAR")
    println("  Optional tablevariable flags for direct gradients (raw mode):")
    println("    --tv-dudz=TABLE.VAR --tv-dthetadz=TABLE.VAR")
    println("  Additional required tablevariable flags for two-level gradients:")
    println("    --tv-u1=TABLE.VAR --tv-u2=TABLE.VAR --tv-theta1=TABLE.VAR --tv-theta2=TABLE.VAR")
end

function parse_flag(args::Vector{String}, key::String, default::String)
    for a in args
        startswith(a, key) || continue
        parts = split(a, "=", limit=2)
        length(parts) == 2 || return default
        return parts[2]
    end
    return default
end

function parse_float_flag(args::Vector{String}, key::String)
    for a in args
        startswith(a, key) || continue
        parts = split(a, "=", limit=2)
        length(parts) == 2 || error("Missing value for $(key)")
        v = tryparse(Float64, parts[2])
        v === nothing && error("Could not parse $(key) value as Float64: $(parts[2])")
        return v
    end
    return nothing
end

function parse_required_flag(args::Vector{String}, key::String)
    value = parse_flag(args, key, "")
    isempty(value) && error("Missing required flag: $(key)=...")
    return value
end

has_flag(args::Vector{String}, flag::String) = any(x -> x == flag, args)

function pick_col(df::DataFrame, aliases::Vector{Symbol}; required::Bool=true)
    colmap = Dict(Symbol(lowercase(String(c))) => c for c in names(df))
    for a in aliases
        key = Symbol(lowercase(String(a)))
        if haskey(colmap, key)
            return colmap[key]
        end
    end
    if required
        error("Missing required column. Tried aliases: $(join(string.(aliases), ", "))")
    end
    return nothing
end

to_float(x) = x isa Missing ? NaN : (x isa Number ? Float64(x) : something(tryparse(Float64, strip(String(x))), NaN))

function col_float(df::DataFrame, col)
    src = df[!, col]
    out = Vector{Float64}(undef, length(src))
    @inbounds for i in eachindex(src)
        out[i] = to_float(src[i])
    end
    return out
end

function maybe_col_float(df::DataFrame, aliases::Vector{Symbol})
    c = pick_col(df, aliases; required=false)
    return c === nothing ? nothing : col_float(df, c)
end

function output_base(output_csv::String)
    return endswith(lowercase(output_csv), ".csv") ? output_csv[1:(end - 4)] : output_csv
end

function write_preprocess_summary(output_csv::String, out::DataFrame; input_source::String, mode::String, phi_target::String, stable_only::Bool, z_m::Float64, d_m::Float64, z_eff::Float64, wtv_eps::Float64, total_rows::Int, neutral_count::Int)
    base = output_base(output_csv)

    rows_written = nrow(out)
    zeta_vals = rows_written > 0 ? Vector{Float64}(out.zeta) : Float64[]
    phi_vals = rows_written > 0 ? Vector{Float64}(out.phi_obs) : Float64[]

    zeta_min = rows_written > 0 ? minimum(zeta_vals) : NaN
    zeta_max = rows_written > 0 ? maximum(zeta_vals) : NaN
    zeta_q05 = rows_written > 0 ? quantile(zeta_vals, 0.05) : NaN
    zeta_q50 = rows_written > 0 ? quantile(zeta_vals, 0.50) : NaN
    zeta_q95 = rows_written > 0 ? quantile(zeta_vals, 0.95) : NaN

    phi_min = rows_written > 0 ? minimum(phi_vals) : NaN
    phi_max = rows_written > 0 ? maximum(phi_vals) : NaN
    phi_q05 = rows_written > 0 ? quantile(phi_vals, 0.05) : NaN
    phi_q50 = rows_written > 0 ? quantile(phi_vals, 0.50) : NaN
    phi_q95 = rows_written > 0 ? quantile(phi_vals, 0.95) : NaN

    stats = DataFrame(
        input_source=[input_source],
        mode=[mode],
        phi_target=[phi_target],
        stable_only=[stable_only],
        z_m=[z_m],
        d_m=[d_m],
        z_eff=[z_eff],
        wtv_eps=[wtv_eps],
        total_rows=[total_rows],
        rows_written=[rows_written],
        neutral_count=[neutral_count],
        zeta_min=[zeta_min],
        zeta_q05=[zeta_q05],
        zeta_q50=[zeta_q50],
        zeta_q95=[zeta_q95],
        zeta_max=[zeta_max],
        phi_min=[phi_min],
        phi_q05=[phi_q05],
        phi_q50=[phi_q50],
        phi_q95=[phi_q95],
        phi_max=[phi_max],
    )
    CSV.write("$(base)_preprocess_stats.csv", stats)

    lines = [
        "# Preprocess Summary",
        "",
        "## Source",
        "",
        "- input source: $(input_source)",
        "- mode: $(mode)",
        "- phi target: $(phi_target)",
        "- stable_only: $(stable_only)",
        "",
        "## Geometry and Thresholds",
        "",
        "- z_m: $(z_m)",
        "- d_m: $(d_m)",
        "- z_eff: $(z_eff)",
        "- wtv_eps: $(wtv_eps)",
        "",
        "## Row Accounting",
        "",
        "- total rows read: $(total_rows)",
        "- rows written: $(rows_written)",
        "- neutral-transition rows flagged: $(neutral_count)",
        "",
        "## Output Distribution",
        "",
        "- zeta min/max: $(zeta_min), $(zeta_max)",
        "- zeta quantiles 5/50/95%: $(zeta_q05), $(zeta_q50), $(zeta_q95)",
        "- phi_obs min/max: $(phi_min), $(phi_max)",
        "- phi_obs quantiles 5/50/95%: $(phi_q05), $(phi_q50), $(phi_q95)",
        "",
        "## Files",
        "",
        "- data: $(output_csv)",
        "- stats csv: $(base)_preprocess_stats.csv",
    ]
    write("$(base)_preprocess_summary.md", join(lines, "\n") * "\n")
end

function get_tablevariable_flags(extra::Vector{String}, mode::String)
    tv = Dict{Symbol, String}()

    tv_mo = parse_flag(extra, "--tv-mo-length", "")
    tv_ustar = parse_flag(extra, "--tv-ustar", "")
    has_direct_mo = !isempty(tv_mo)

    if has_direct_mo
        tv[:mo_length] = tv_mo
        if !isempty(tv_ustar)
            tv[:u_star] = tv_ustar
        end
    else
        tv[:uw] = parse_required_flag(extra, "--tv-uw")
        tv[:vw] = parse_required_flag(extra, "--tv-vw")
        tv[:wthetav] = parse_required_flag(extra, "--tv-wthetav")
        tv[:thetav] = parse_required_flag(extra, "--tv-thetav")
    end

    tv_dudz = parse_flag(extra, "--tv-dudz", "")
    tv_dthetadz = parse_flag(extra, "--tv-dthetadz", "")
    if !isempty(tv_dudz)
        tv[:dudz] = tv_dudz
    end
    if !isempty(tv_dthetadz)
        tv[:dthetadz] = tv_dthetadz
    end

    if mode == "two-level"
        tv[:u1] = parse_required_flag(extra, "--tv-u1")
        tv[:u2] = parse_required_flag(extra, "--tv-u2")
        tv[:theta1] = parse_required_flag(extra, "--tv-theta1")
        tv[:theta2] = parse_required_flag(extra, "--tv-theta2")
    end

    return tv
end

function fetch_smear_dataframe(station_label::String, extra::Vector{String}, mode::String)
    from_iso = parse_required_flag(extra, "--from")
    to_iso = parse_required_flag(extra, "--to")
    interval = parse_flag(extra, "--interval", "30")
    aggregation = parse_flag(extra, "--aggregation", "ARITHMETIC")
    quality = parse_flag(extra, "--quality", "ANY")

    tv = get_tablevariable_flags(extra, mode)

    tablevars = String[]
    for k in keys(tv)
        push!(tablevars, tv[k])
    end
    tablevars = unique(tablevars)

    params = String[]
    push!(params, "from=$(from_iso)")
    push!(params, "to=$(to_iso)")
    push!(params, "interval=$(interval)")
    push!(params, "aggregation=$(aggregation)")
    push!(params, "quality=$(quality)")
    for tv_name in tablevars
        push!(params, "tablevariable=$(tv_name)")
    end

    url = "https://smear-backend-avaa-smear-prod.2.rahtiapp.fi/search/timeseries/csv?" * join(params, "&")
    tmp = Downloads.download(url)
    df_api = CSV.read(tmp, DataFrame)
    rm(tmp; force=true)

    has_stamp = all(col -> col in names(df_api), [:Year, :Month, :Day, :Hour, :Minute, :Second])
    if has_stamp
        ts = Vector{String}(undef, nrow(df_api))
        for i in 1:nrow(df_api)
            y = Int(to_float(df_api.Year[i]))
            mo = Int(to_float(df_api.Month[i]))
            d = Int(to_float(df_api.Day[i]))
            h = Int(to_float(df_api.Hour[i]))
            mi = Int(to_float(df_api.Minute[i]))
            s = Int(to_float(df_api.Second[i]))
            dt = DateTime(y, mo, d, h, mi, s)
            ts[i] = Dates.format(dt, dateformat"yyyy-mm-ddTHH:MM:SS")
        end
        df_api.timestamp = ts
    end

    colmap = Dict(lowercase(String(c)) => c for c in names(df_api))
    for (k, tv_name) in tv
        key = lowercase(tv_name)
        haskey(colmap, key) || error("Requested tablevariable $(tv_name) not found in API response.")
        col = colmap[key]
        df_api[!, k] = df_api[!, col]
    end

    println("Fetched $(nrow(df_api)) rows from SmartSMEAR for station label $(station_label)")
    return df_api
end

function main()
    if length(ARGS) < 4
        print_usage()
        return
    end

    input_csv = ARGS[1]
    output_csv = ARGS[2]
    z_m = parse(Float64, ARGS[3])
    d_m = parse(Float64, ARGS[4])

    extra = length(ARGS) > 4 ? ARGS[5:end] : String[]
    stable_only = has_flag(extra, "--stable-only")
    phi_target = String(parse_flag(extra, "--phi", "phi_m"))
    wtv_eps = parse_float_flag(extra, "--wtv-eps")
    wtv_eps = wtv_eps === nothing ? DEFAULT_WTV_EPS : wtv_eps
    mode = String(parse_flag(extra, "--mode", "raw"))
    if !(phi_target in ("phi_m", "phi_h"))
        error("--phi must be phi_m or phi_h")
    end
    if !(mode in ("raw", "two-level", "api-smear"))
        error("--mode must be raw, two-level, or api-smear")
    end

    z_eff = z_m - d_m
    z_eff > 0 || error("Need z_m - d_m > 0. Got z_m=$(z_m), d_m=$(d_m)")

    df = if mode == "api-smear"
        station_label = input_csv
        profile_mode = String(parse_flag(extra, "--profile-mode", "raw"))
        if !(profile_mode in ("raw", "two-level"))
            error("For --mode=api-smear, --profile-mode must be raw or two-level")
        end
        mode = profile_mode
        fetch_smear_dataframe(station_label, extra, mode)
    else
        CSV.read(input_csv, DataFrame)
    end

    c_mo = pick_col(df, [:mo_length, :MO_length, :L, :monin_obukhov_length]; required=false)
    c_ustar_direct = pick_col(df, [:u_star, :ustar, :uStar]; required=false)
    have_direct_mo = c_mo !== nothing

    c_uw = pick_col(df, [:uw, :u_w_cov, :u_prime_w_prime]; required=!have_direct_mo)
    c_vw = pick_col(df, [:vw, :v_w_cov, :v_prime_w_prime]; required=!have_direct_mo)
    c_wtv = pick_col(df, [:wthetav, :w_theta_v, :wthetav_cov]; required=!have_direct_mo)
    c_thetav = pick_col(df, [:thetav, :theta_v, :theta_v_mean]; required=!have_direct_mo)

    uw = c_uw === nothing ? fill(NaN, nrow(df)) : col_float(df, c_uw)
    vw = c_vw === nothing ? fill(NaN, nrow(df)) : col_float(df, c_vw)
    wtv = c_wtv === nothing ? fill(NaN, nrow(df)) : col_float(df, c_wtv)
    thetav = c_thetav === nothing ? fill(NaN, nrow(df)) : col_float(df, c_thetav)
    mo_length = c_mo === nothing ? fill(NaN, nrow(df)) : col_float(df, c_mo)
    ustar_direct = c_ustar_direct === nothing ? fill(NaN, nrow(df)) : col_float(df, c_ustar_direct)

    dudz = nothing
    dthetadz = nothing

    if mode == "raw"
        dudz = maybe_col_float(df, [:dudz, :dU_dz, :dUdz])
        dthetadz = maybe_col_float(df, [:dthetadz, :dtheta_dz, :dthetav_dz])
    else
        z1 = parse_float_flag(extra, "--z1")
        z2 = parse_float_flag(extra, "--z2")
        (z1 === nothing || z2 === nothing) && error("Mode two-level requires --z1 and --z2 flags.")
        dz = z2 - z1
        abs(dz) > 0 || error("Mode two-level requires z2 != z1.")

        c_u1 = pick_col(df, [:u1, :u_low, :U1, :wind_low])
        c_u2 = pick_col(df, [:u2, :u_high, :U2, :wind_high])
        c_t1 = pick_col(df, [:theta1, :theta_low, :T1, :temp_low])
        c_t2 = pick_col(df, [:theta2, :theta_high, :T2, :temp_high])

        u1 = col_float(df, c_u1)
        u2 = col_float(df, c_u2)
        t1 = col_float(df, c_t1)
        t2 = col_float(df, c_t2)

        dudz = Vector{Float64}(undef, length(u1))
        dthetadz = Vector{Float64}(undef, length(t1))
        @inbounds for i in eachindex(u1)
            dudz[i] = (u2[i] - u1[i]) / dz
            dthetadz[i] = (t2[i] - t1[i]) / dz
        end
    end

    n = nrow(df)
    ustar = similar(uw)
    L = similar(uw)
    zeta = similar(uw)
    neutral_flag = falses(n)
    phi_m = fill(NaN, n)
    phi_h = fill(NaN, n)

    for i in 1:n
        tau = sqrt(uw[i]^2 + vw[i]^2)
        ustar_from_tau = tau^(0.5)
        if isfinite(ustar_direct[i]) && ustar_direct[i] > 0
            ustar[i] = ustar_direct[i]
        else
            ustar[i] = ustar_from_tau
        end

        if have_direct_mo && isfinite(mo_length[i]) && mo_length[i] != 0.0
            L[i] = mo_length[i]
            zeta[i] = z_eff / L[i]
        else
            if !(isfinite(ustar[i]) && isfinite(thetav[i]) && isfinite(wtv[i]) && ustar[i] > 0 && thetav[i] > 0)
                L[i] = NaN
                zeta[i] = NaN
                continue
            end

            if abs(wtv[i]) <= wtv_eps
                # Treat near-zero buoyancy flux as neutral transition to avoid dropping rows.
                neutral_flag[i] = true
                zeta[i] = 0.0
                L[i] = wtv[i] < 0 ? Inf : -Inf
            else
                L[i] = -(ustar[i]^3 * thetav[i]) / (KAPPA * G * wtv[i])
                zeta[i] = z_eff / L[i]
            end
        end

        if dudz !== nothing && isfinite(dudz[i]) && ustar[i] > 0
            phi_m[i] = KAPPA * z_eff * dudz[i] / ustar[i]
        end

        if dthetadz !== nothing && isfinite(dthetadz[i])
            theta_star = -wtv[i] / ustar[i]
            if isfinite(theta_star) && abs(theta_star) > 1e-12
                phi_h[i] = KAPPA * z_eff * dthetadz[i] / theta_star
            end
        end
    end

    quality_pass = isfinite.(zeta) .& isfinite.(ustar)
    if stable_only
        quality_pass .&= zeta .> 0
    end

    phi_obs = phi_target == "phi_m" ? phi_m : phi_h
    quality_pass .&= isfinite.(phi_obs)

    out = DataFrame()
    if :time in names(df)
        out.time = df.time
    elseif :timestamp in names(df)
        out.time = df.timestamp
    else
        out.time = collect(1:n)
    end

    out.zeta = zeta
    out.phi_obs = phi_obs
    out.phi_m = phi_m
    out.phi_h = phi_h
    out.u_star = ustar
    out.L = L
    out.neutral_transition = neutral_flag
    out.quality_pass = quality_pass

    out = out[out.quality_pass .== true, :]
    CSV.write(output_csv, out)
    write_preprocess_summary(
        output_csv,
        out;
        input_source=input_csv,
        mode=mode,
        phi_target=phi_target,
        stable_only=stable_only,
        z_m=z_m,
        d_m=d_m,
        z_eff=z_eff,
        wtv_eps=wtv_eps,
        total_rows=n,
        neutral_count=count(neutral_flag),
    )

    println("Wrote $(nrow(out)) filtered rows to $(output_csv)")
    println("phi_obs source: $(phi_target)")
    println("stable_only: $(stable_only)")
    println("mode: $(mode)")
    println("wtv_eps: $(wtv_eps)")
end

try
    main()
catch err
    println(stderr, "Error: ", err)
    exit(1)
end
