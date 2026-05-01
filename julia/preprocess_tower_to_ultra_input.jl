#!/usr/bin/env julia

using CSV, DataFrames

const KAPPA = 0.4
const G = 9.81

function print_usage()
    println("Usage:")
    println("  julia julia/preprocess_tower_to_ultra_input.jl <input_csv> <output_csv> <z_m> <d_m> [--stable-only] [--phi=phi_m|phi_h] [--mode=raw|two-level] [--z1=<m> --z2=<m>]")
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

function to_float(x)
    if x isa Missing
        return NaN
    elseif x isa Number
        return Float64(x)
    elseif x isa AbstractString
        y = tryparse(Float64, strip(x))
        return y === nothing ? NaN : y
    else
        return NaN
    end
end

function col_float(df::DataFrame, col::Symbol)
    return [to_float(v) for v in df[!, col]]
end

function maybe_col_float(df::DataFrame, aliases::Vector{Symbol})
    c = pick_col(df, aliases; required=false)
    return c === nothing ? nothing : col_float(df, c)
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
    phi_target = parse_flag(extra, "--phi", "phi_m")
    mode = parse_flag(extra, "--mode", "raw")
    if !(phi_target in ("phi_m", "phi_h"))
        error("--phi must be phi_m or phi_h")
    end
    if !(mode in ("raw", "two-level"))
        error("--mode must be raw or two-level")
    end

    z_eff = z_m - d_m
    z_eff > 0 || error("Need z_m - d_m > 0. Got z_m=$(z_m), d_m=$(d_m)")

    df = CSV.read(input_csv, DataFrame)

    c_uw = pick_col(df, [:uw, :u_w_cov, :u_prime_w_prime])
    c_vw = pick_col(df, [:vw, :v_w_cov, :v_prime_w_prime])
    c_wtv = pick_col(df, [:wthetav, :w_theta_v, :wthetav_cov])
    c_thetav = pick_col(df, [:thetav, :theta_v, :theta_v_mean])

    uw = col_float(df, c_uw)
    vw = col_float(df, c_vw)
    wtv = col_float(df, c_wtv)
    thetav = col_float(df, c_thetav)

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

        dudz = (u2 .- u1) ./ dz
        dthetadz = (t2 .- t1) ./ dz
    end

    n = nrow(df)
    ustar = similar(uw)
    L = similar(uw)
    zeta = similar(uw)
    phi_m = fill(NaN, n)
    phi_h = fill(NaN, n)

    for i in 1:n
        tau = sqrt(uw[i]^2 + vw[i]^2)
        ustar[i] = tau^(0.5)

        if !(isfinite(ustar[i]) && isfinite(thetav[i]) && isfinite(wtv[i]) && ustar[i] > 0 && thetav[i] > 0 && abs(wtv[i]) > 1e-12)
            L[i] = NaN
            zeta[i] = NaN
            continue
        end

        L[i] = -(ustar[i]^3 * thetav[i]) / (KAPPA * G * wtv[i])
        zeta[i] = z_eff / L[i]

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
    out.quality_pass = quality_pass

    out = out[out.quality_pass .== true, :]
    CSV.write(output_csv, out)

    println("Wrote $(nrow(out)) filtered rows to $(output_csv)")
    println("phi_obs source: $(phi_target)")
    println("stable_only: $(stable_only)")
    println("mode: $(mode)")
end

main()
