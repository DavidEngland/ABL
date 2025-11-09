module MOSTProfiles

export make_profile,
       F_from, ri_from_zeta,
       zeta_from_ri_series, zeta_from_ri_newton,
       ri_to_phi_wrappers,
       compute_vlog_wlog, rig_derivatives_zeta,
       dzeta_terms, map_curvature_z, omit_error_metric, map_curvature_auto, map_curvature_point

# --------------------------- Utilities -----------------------------------------

const TFloat = Union{Float32, Float64}
const RealOrVec = Union{Real, AbstractVector{<:Real}}

@inline function centraldiff(f::Function, x::Real; h::Real=1e-6)
    return (f(x + h) - f(x - h)) / (2h)
end

@inline function seconddiff(f::Function, x::Real; h::Real=1e-6)
    return (f(x + h) - 2f(x) + f(x - h)) / (h * h)
end

# Vectorized wrapper
centraldiff(f::Function, x::AbstractVector; h::Real=1e-6) = [centraldiff(f, xi; h=h) for xi in x]
seconddiff(f::Function, x::AbstractVector; h::Real=1e-6)  = [seconddiff(f, xi; h=h) for xi in x]

# Safe division
@inline safediv(a::Real, b::Real; eps::Real=1e-300) = a / (abs(b) > eps ? b : eps)

# --------------------------- Profile registry ----------------------------------

"""
    make_profile(tag::AbstractString, pars::Dict)

Return (phi_m(ζ), phi_h(ζ)) callables for the requested profile family.

Supported tags:
- "BD_PL"      : Power-law Businger–Dyer (am,bm,ah,bh)
- "BD_CLASSIC" : Classical composite (a,cm,ch,pm_exp,ph_exp)
- "HOG88"      : Linear stable (cm,ch,c0h)
- "QSBL"       : Quadratic stable surrogate (am,bm,ah,bh) for ζ≥0
- "CB"         : Cheng–Brutsaert monotone (gm,pm,gh,ph)
- "RPL"        : Regularized power law (alpha_m,beta_m,delta_m,alpha_h,beta_h,delta_h)
- "VEXP"       : Variable exponent (alpha_m,beta_m,eta_m,alpha_h,beta_h,eta_h)
- "DTP"        : Dynamic turbulent Prandtl wrapper; returns base φ (base_tag, base_pars)
- "URC"        : Ri-based closure; returns functions of Ri (b_m,Ri_c,e_m[,b_h,e_h,Ri_c_h])

Notes
- Domain guard: power-like forms require 1 - βζ > eps.
- For "URC", functions expect Ri (handled via `ri_to_phi_wrappers`).
"""
function make_profile(tag::AbstractString, pars::Dict)
    t = uppercase(tag)

    if t == "BD_PL"
        am,bm,ah,bh = pars[:am], pars[:bm], pars[:ah], pars[:bh]
        ϵ = get(pars, :eps, 1e-12)
        ϕm = z -> begin
            d = 1 - bm*z
            d > ϵ ? d^(-am) : NaN
        end
        ϕh = z -> begin
            d = 1 - bh*z
            d > ϵ ? d^(-ah) : NaN
        end
        return ϕm, ϕh
    end

    if t == "BD_CLASSIC"
        a   = get(pars, :a, 16.0)
        cm  = get(pars, :cm, 5.0)
        ch  = get(pars, :ch, 7.0)
        pme = get(pars, :pm_exp, -0.25)
        phe = get(pars, :ph_exp, -0.5)
        ϕm = z -> z < 0 ? (1 - a*z)^(pme) : 1 + cm*z
        ϕh = z -> z < 0 ? (1 - a*z)^(phe) : 1 + ch*z
        return ϕm, ϕh
    end

    if t == "HOG88"
        cm  = get(pars, :cm, 5.0)
        ch  = get(pars, :ch, 7.8)
        c0h = get(pars, :c0h, 0.95)
        return z -> 1 + cm*z, z -> c0h + ch*z
    end

    if t == "QSBL"
        am,bm,ah,bh = pars[:am], pars[:bm], pars[:ah], pars[:bh]
        ϕm = z -> 1 + am*z + bm*z*z
        ϕh = z -> 1 + ah*z + bh*z*z
        return ϕm, ϕh
    end

    if t == "CB"
        gm,pm,gh,ph = pars[:gm], pars[:pm], pars[:gh], pars[:ph]
        return z -> (1 + gm*abs(z))^pm,
               z -> (1 + gh*abs(z))^ph
    end

    if t == "RPL"
        am,bm,dm = pars[:alpha_m], pars[:beta_m], pars[:delta_m]
        ah,bh,dh = pars[:alpha_h], pars[:beta_h], pars[:delta_h]
        g(b,d,z) = (b*z) / (1 + d*b*z)
        return z -> (1 + g(bm,dm,z))^am,
               z -> (1 + g(bh,dh,z))^ah
    end

    if t == "VEXP"
        am,bm,em = pars[:alpha_m], pars[:beta_m], pars[:eta_m]
        ah,bh,eh = pars[:alpha_h], pars[:beta_h], pars[:eta_h]
        ϵ = get(pars, :eps, 1e-12)
        ϕm = z -> begin
            d = 1 - bm*z
            d > ϵ ? d^(-am*(1 + em*z)) : NaN
        end
        ϕh = z -> begin
            d = 1 - bh*z
            d > ϵ ? d^(-ah*(1 + eh*z)) : NaN
        end
        return ϕm, ϕh
    end

    if t == "DTP"
        base_tag  = pars[:base_tag]
        base_pars = pars[:base_pars]
        return make_profile(base_tag, base_pars)
    end

    if t == "URC"
        b_m, Ri_c, e_m = pars[:b_m], pars[:Ri_c], pars[:e_m]
        fm = Ri -> (1 + b_m * Ri / Ri_c)^(-e_m)
        fh = haskey(pars, :b_h) && haskey(pars, :e_h) ? (Ri -> (1 + pars[:b_h] * Ri / get(pars, :Ri_c_h, Ri_c))^(-pars[:e_h])) : (Ri -> NaN)
        return fm, fh
    end

    error("unknown profile tag '$tag'")
end

# --------------------------- Ri / ζ utilities ----------------------------------

F_from(ϕm::Function, ϕh::Function) = (z -> ϕh(z) / (ϕm(z)^2))

ri_from_zeta(ζ::Real, ϕm::Function, ϕh::Function) = ζ * F_from(ϕm, ϕh)(ζ)
ri_from_zeta(ζ::AbstractVector, ϕm::Function, ϕh::Function) = ζ .* (F_from(ϕm, ϕh)).(ζ)

"""
    zeta_from_ri_series(Ri, Δ, c1)

Near-neutral inversion ζ(Ri) ≈ Ri - Δ Ri² + (1.5Δ² - 0.5 c1) Ri³.
"""
zeta_from_ri_series(Ri::Real, Δ::Real, c1::Real) = Ri - Δ*Ri*Ri + (1.5*Δ*Δ - 0.5*c1)*(Ri^3)

"""
    zeta_from_ri_newton(Ri_target, ϕm, ϕh, z0; tol=1e-10, maxit=20)

Refine ζ solving f(ζ)=ζ F(ζ) - Ri_target = 0 via Newton. Uses V_log = d(log F)/dζ.
"""
function zeta_from_ri_newton(Ri_target::Real, ϕm::Function, ϕh::Function, z0::Real; tol::Real=1e-10, maxit::Int=20)
    F = F_from(ϕm, ϕh)
    z = z0
    @inbounds for _ in 1:maxit
        Vlog = centraldiff(zz -> log(F(zz)), z)
        f  = z*F(z) - Ri_target
        Fz = F(z)
        fp = Fz + z*Fz*Vlog  # F(1 + z V_log)
        if fp == 0.0
            break
        end
        dz = f / fp
        z -= dz
        if abs(dz) < tol
            break
        end
    end
    return z
end

"""
    compute_vlog_wlog(ζ, ϕm, ϕh; h=1e-6)

Return (v_m, v_h, V_log, W_log) at ζ using central differences.
"""
function compute_vlog_wlog(ζ::Real, ϕm::Function, ϕh::Function; h::Real=1e-6)
    pm = ϕm(ζ); ph = ϕh(ζ)
    vm = safediv(centraldiff(ϕm, ζ; h=h), pm)
    vh = safediv(centraldiff(ϕh, ζ; h=h), ph)
    V_log = vh - 2vm
    # Differentiate V_log directly for W_log
    V_of = zz -> begin
        pmt = ϕm(zz); pht = ϕh(zz)
        vmₜ = safediv(centraldiff(ϕm, zz; h=h), pmt)
        vhₜ = safediv(centraldiff(ϕh, zz; h=h), pht)
        vhₜ - 2vmₜ
    end
    W_log = centraldiff(V_of, ζ; h=h)
    return vm, vh, V_log, W_log
end

"""
    rig_derivatives_zeta(ζ, ϕm, ϕh)

Return (dRi/dζ, d²Ri/dζ², F, V_log) at ζ using unified log-derivatives:
  F = φ_h / φ_m²
  dRi/dζ = F (1 + ζ V_log)
  d²Ri/dζ² = F [ 2 V_log + ζ(V_log² - W_log) ]
"""
function rig_derivatives_zeta(ζ::Real, ϕm::Function, ϕh::Function)
    F = F_from(ϕm, ϕh)(ζ)
    _, _, Vlog, Wlog = compute_vlog_wlog(ζ, ϕm, ϕh)
    dRi_dζ  = F * (1 + ζ * Vlog)
    d2Ri_dζ2 = F * (2Vlog + ζ * (Vlog*Vlog - Wlog))
    return dRi_dζ, d2Ri_dζ2, F, Vlog
end

# Vectorized variant
function rig_derivatives_zeta(ζ::AbstractVector, ϕm::Function, ϕh::Function)
    d1 = similar(ζ); d2 = similar(ζ); Fv = similar(ζ); Vv = similar(ζ)
    @inbounds for i in eachindex(ζ)
        d1[i], d2[i], Fv[i], Vv[i] = rig_derivatives_zeta(ζ[i], ϕm, ϕh)
    end
    return d1, d2, Fv, Vv
end

"""
    ri_to_phi_wrappers(tag, pars; Δ=nothing, c1=nothing)

Return (f_m(Ri), f_h(Ri)) closures from ζ-space profile or pass-through for URC.
For DTP, heat closure uses Pr_t(Ri) = 1 + a1 Ri + a2 Ri² applied to base φ_m.
"""
function ri_to_phi_wrappers(tag::AbstractString, pars::Dict; Δ=nothing, c1=nothing)
    if uppercase(tag) == "URC"
        return make_profile(tag, pars)
    end
    ϕm, ϕh = make_profile(tag, pars)

    zeta_of_Ri = Ri -> begin
        z0 = (Δ !== nothing && c1 !== nothing) ? zeta_from_ri_series(Ri, Δ, c1) : Ri
        zeta_from_ri_newton(Ri, ϕm, ϕh, z0)
    end

    if uppercase(tag) == "DTP"
        base_tag  = pars[:base_tag]
        base_pars = pars[:base_pars]
        a1 = get(pars, :a1, 0.0)
        a2 = get(pars, :a2, 0.0)
        base_m, _ = make_profile(base_tag, base_pars)
        fm = Ri -> begin
            ζ = zeta_of_Ri(Ri)
            base_m(ζ)
        end
        fh = Ri -> begin
            ζ = zeta_of_Ri(Ri)
            Prt = 1 + a1*Ri + a2*Ri*Ri
            Prt * base_m(ζ)
        end
        return fm, fh
    end

    fm = Ri -> begin
        ζ = zeta_of_Ri(Ri)
        ϕm(ζ)
    end
    fh = Ri -> begin
        ζ = zeta_of_Ri(Ri)
        ϕh(ζ)
    end
    return fm, fh
end

# --------------------------- Variable-L height mapping --------------------------

"""
    dzeta_terms(z, L)

Compute dζ/dz and d²ζ/dz² for ζ = z / L(z) given arrays z, L(z).
"""
function dzeta_terms(z::AbstractVector, L::AbstractVector)
    @assert length(z) == length(L)
    L1 = gradient(L, z)
    L2 = gradient(L1, z)
    dzeta_dz   = (L .- z .* L1) ./ (L .* L)
    d2zeta_dz2 = (-2 .* L1) ./ (L .* L) .- (z .* L2) ./ (L .* L) .+ 2 .* z .* (L1 .* L1) ./ (L .* L .* L)
    return dzeta_dz, d2zeta_dz2
end

# Simple centered gradient on irregular grid
function gradient(y::AbstractVector, x::AbstractVector)
    n = length(y)
    g = similar(y)
    @inbounds begin
        if n >= 3
            # endpoints: second-order one-sided
            g[1] = ( -3y[1] + 4y[2] - y[3]) / (x[3] - x[1])
            g[end] = ( 3y[end] - 4y[end-1] + y[end-2]) / (x[end] - x[end-2])
            for i in 2:n-1
                g[i] = (y[i+1] - y[i-1]) / (x[i+1] - x[i-1])
            end
        elseif n == 2
            g[1] = (y[2] - y[1]) / (x[2] - x[1])
            g[2] = g[1]
        else
            g[1] = 0
        end
    end
    return g
end

"""
    map_curvature_z(z, L, dRi_dζ, d2Ri_dζ2)

Map ζ-space derivatives to height-space curvature:
  ∂²Ri/∂z² = (dζ/dz)² ∂²Ri/∂ζ² + (d²ζ/dz²) ∂Ri/∂ζ
"""
function map_curvature_z(z::AbstractVector, L::AbstractVector,
                         dRi_dζ::AbstractVector, d2Ri_dζ2::AbstractVector)
    dzeta_dz, d2zeta_dz2 = dzeta_terms(z, L)
    return (dzeta_dz .* dzeta_dz) .* d2Ri_dζ2 .+ d2zeta_dz2 .* dRi_dζ
end

"""
    omit_error_metric(dζdz, d²ζdz², dRi_dζ, d²Ri_dζ²)

E_omit = |(d²ζ/dz² * dRi/dζ) / ((dζ/dz)² * d²Ri/dζ²)|
If E_omit < eps, constant-L shortcut is acceptable at that level.
"""
function omit_error_metric(dζdz::AbstractVector, d²ζdz²::AbstractVector,
                           dRi_dζ::AbstractVector, d²Ri_dζ²::AbstractVector)
    denom = (dζdz .* dζdz) .* d²Ri_dζ²
    safe = map(x -> abs(x) > 0 ? x : NaN, denom)
    return abs.(d²ζdz² .* dRi_dζ) ./ abs.(safe)
end

"""
    map_curvature_auto(z, L, dRi_dζ, d²Ri_dζ²; eps=0.05)

Auto-select constant-L or full mapping per level using E_omit threshold.
"""
function map_curvature_auto(z::AbstractVector, L::AbstractVector,
                            dRi_dζ::AbstractVector, d²Ri_dζ²::AbstractVector; eps::Real=0.05)
    dζdz, d²ζdz² = dzeta_terms(z, L)
    E = omit_error_metric(dζdz, d²ζdz², dRi_dζ, d²Ri_dζ²)
    curv_const = d²Ri_dζ² ./ (L .* L)
    curv_full  = (dζdz .* dζdz) .* d²Ri_dζ² .+ d²ζdz² .* dRi_dζ
    use_const = map(e -> !isfinite(e) || e < eps, E)
    return ifelse.(use_const, curv_const, curv_full)
end

"""
    map_curvature_point(z, L, dL, d2L, dRi_dζ, d²Ri_dζ²)

Pointwise mapping when L, dL, d²L are known at z.
"""
function map_curvature_point(z::Real, L::Real, dL::Real, d2L::Real,
                             dRi_dζ::Real, d²Ri_dζ²::Real)
    dζdz   = (L - z*dL) / (L * L)
    d2ζdz2 = (-2dL) / (L * L) - (z*d2L) / (L * L) + 2z*(dL*dL) / (L * L * L)
    return (dζdz*dζdz) * d²Ri_dζ² + d2ζdz2 * dRi_dζ
end

# --------------------------- Minimal example (comment) --------------------------
# Example:
# using .MOSTProfiles
# pars = Dict(:am=>0.5,:bm=>16.0,:ah=>0.5,:bh=>16.0)
# ϕm, ϕh = make_profile("BD_PL", pars)
# ζ = 0.03
# d1, d2, F, V = rig_derivatives_zeta(ζ, ϕm, ϕh)
# Ri = ri_from_zeta(ζ, ϕm, ϕh)
# ζ0 = zeta_from_ri_series(Ri, pars[:ah]*pars[:bh]-2*pars[:am]*pars[:bm],
#                               pars[:ah]*pars[:bh]^2-2*pars[:am]*pars[:bm]^2)
# ζ_refined = zeta_from_ri_newton(Ri, ϕm, ϕh, ζ0)

end # module
