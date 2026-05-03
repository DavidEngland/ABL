#!/usr/bin/env julia
# preprocess_sheba_main.jl
# ------------------------------------------------------------------
# Download (or read) the NOAA PSL SHEBA ASFG main_file6_hd.txt and
# compute zeta = z/L and phi_m from the 2-level wind profile.
#
# Output: CSV with columns  time, zeta, phi_obs
#
# Usage:
#   julia julia/preprocess_sheba_main.jl <output_csv> [input_txt]
#
# If input_txt is omitted the file is downloaded from NOAA PSL.
# ------------------------------------------------------------------

using CSV, DataFrames, Statistics, Downloads

const DATA_URL = "https://psl.noaa.gov/arctic/sheba/netcdf/PerssonDatasets/main_file6_hd.txt"
const KAPPA    = 0.40     # von Kármán constant
const G        = 9.81     # m s⁻²
const RHO_CP   = 1200.0   # J m⁻³ K⁻¹  (ρ cₚ for arctic surface air)
const MISSING_VAL = 9000.0  # values ≥ 9000 are flagged missing (9999)

# Reference heights (metres) for the 2-level wind profile
const Z_LO = 2.5
const Z_HI = 10.0
const Z_REF = sqrt(Z_LO * Z_HI)   # geometric-mean reference height ≈ 5.0 m

# Quality-control thresholds
const USTAR_MIN   = 0.05   # m s⁻¹
const HS_ABS_MIN  = 2.0    # W m⁻²  (avoid near-neutral scatter)
const ZETA_MAX    = 10.0   # upper limit for stable zeta
const PHI_MAX     = 30.0   # upper limit on dimensionless wind shear

function compute_L(ustar, hs_wm2, T_C)
    # Obukhov length [m]  (positive in stable stratification)
    # L = - u*³ Tᵥ / (κ g H)   where H = hs/(ρ cₚ) [K m s⁻¹]
    T_K = T_C + 273.15
    H   = hs_wm2 / RHO_CP           # kinematic heat flux [K m s⁻¹]
    if abs(H) < 1e-10 || ustar < 1e-4
        return NaN
    end
    return -(ustar^3 * T_K) / (KAPPA * G * H)
end

function compute_phi_m(ws_lo, ws_hi, ustar)
    # Dimensionless wind shear at geometric-mean height Z_REF
    dU_dz = (ws_hi - ws_lo) / (Z_HI - Z_LO)
    return (KAPPA * Z_REF / ustar) * dU_dz
end

function load_sheba(path::AbstractString)
    # File is tab-delimited with a single header row; missing = 9999
    df = CSV.read(path, DataFrame; delim='\t', header=1,
                  missingstring=["9999", "999"],
                  normalizenames=true)
    return df
end

function main()
    if length(ARGS) < 1
        println("Usage: julia preprocess_sheba_main.jl <output_csv> [input_txt]")
        println("  output_csv : path for the output zeta/phi_obs CSV")
        println("  input_txt  : path to local main_file6_hd.txt (downloaded if omitted)")
        exit(1)
    end

    out_csv   = ARGS[1]
    input_txt = length(ARGS) >= 2 ? ARGS[2] : nothing

    # Ensure output directory exists
    mkpath(dirname(out_csv))

    # Download if needed
    local_txt = if input_txt !== nothing
        input_txt
    else
        tmp = tempname() * ".txt"
        @info "Downloading SHEBA main_file6_hd.txt from NOAA PSL …"
        Downloads.download(DATA_URL, tmp)
        @info "Downloaded to $tmp"
        tmp
    end

    @info "Reading $local_txt …"
    df = load_sheba(local_txt)
    @info "Loaded $(nrow(df)) rows, $(ncol(df)) columns"

    # Identify columns robustly (names may be whitespace-trimmed after normalizenames)
    # normalizenames converts e.g. "u*" → "u_" but let's check
    col_names = names(df)
    @info "Columns: $(join(col_names, ", "))"

    # Map expected columns
    function get_col(df, candidates...)
        for c in candidates
            if c in col_names
                return df[!, c]
            end
        end
        error("None of $(candidates) found in columns: $(col_names)")
    end

    jd     = get_col(df, "JD", "JJD")
    ws_lo  = get_col(df, "ws2_5", "ws2.5")
    ws_hi  = get_col(df, "ws10")
    T_lo   = get_col(df, "T2_5", "T2.5")
    T_hi   = get_col(df, "T10")
    ustar  = get_col(df, "u_", "u*", "ustar")
    hs     = get_col(df, "hs")

    n_total = nrow(df)

    # Build output arrays
    times    = Float64[]
    zeta_out = Float64[]
    phi_out  = Float64[]

    n_missing_ustar = 0
    n_missing_ws    = 0
    n_missing_hs    = 0
    n_unstable      = 0
    n_qc_fail       = 0
    n_good          = 0

    for i in 1:n_total
        # Extract values; treat missing or ≥ MISSING_VAL as bad
        u  = coalesce(ustar[i], NaN)
        h  = coalesce(hs[i],    NaN)
        wl = coalesce(ws_lo[i], NaN)
        wh = coalesce(ws_hi[i], NaN)
        tl = coalesce(T_lo[i],  NaN)
        th = coalesce(T_hi[i],  NaN)
        t  = coalesce(jd[i],    NaN)

        if isnan(u) || u >= MISSING_VAL;  n_missing_ustar += 1; continue; end
        if isnan(h) || abs(h) >= MISSING_VAL; n_missing_hs += 1; continue; end
        if isnan(wl) || wl >= MISSING_VAL || isnan(wh) || wh >= MISSING_VAL
            n_missing_ws += 1; continue
        end
        if isnan(tl) || tl >= MISSING_VAL || isnan(th) || th >= MISSING_VAL
            n_missing_ws += 1; continue
        end

        # Basic QC
        if u < USTAR_MIN;                n_qc_fail += 1; continue; end
        if abs(h) < HS_ABS_MIN;          n_qc_fail += 1; continue; end
        if wh < wl;                      n_qc_fail += 1; continue; end  # non-monotonic profile

        T_avg = 0.5 * (tl + th)
        L = compute_L(u, h, T_avg)
        if isnan(L) || isinf(L);         n_qc_fail += 1; continue; end

        zeta = Z_REF / L

        # Stable: zeta > 0
        if zeta <= 0.0;    n_unstable += 1; continue; end
        if zeta > ZETA_MAX; n_qc_fail += 1; continue; end

        phi = compute_phi_m(wl, wh, u)
        if isnan(phi) || isinf(phi) || phi <= 0.0 || phi > PHI_MAX
            n_qc_fail += 1; continue
        end

        push!(times,    isnan(t) ? NaN : t)
        push!(zeta_out, zeta)
        push!(phi_out,  phi)
        n_good += 1
    end

    @info "Preprocessing summary:"
    @info "  Total rows           : $n_total"
    @info "  Missing u*           : $n_missing_ustar"
    @info "  Missing hs           : $n_missing_hs"
    @info "  Missing wind/T       : $n_missing_ws"
    @info "  Unstable (zeta ≤ 0)  : $n_unstable"
    @info "  QC failures          : $n_qc_fail"
    @info "  Good stable rows     : $n_good"

    if n_good < 10
        @warn "Very few good rows ($n_good). Check data format."
    end

    out_df = DataFrame(time=times, zeta=zeta_out, phi_obs=phi_out)
    CSV.write(out_csv, out_df)
    @info "Written to $out_csv"
end

main()
