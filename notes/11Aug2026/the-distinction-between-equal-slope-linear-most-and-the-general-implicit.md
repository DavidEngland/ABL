The distinction between equal-slope linear MOST and the general implicit relationship $\zeta(Ri_g)$ tightens the mathematical foundation, while the ClosureState abstraction provides an extensible software boundary.  
  
## 1. Mathematical Solution for General Slopes ($\beta_m \neq \beta_h$)  
For general linear flux-profile relations $\phi_m(\zeta) = 1 + \beta_m \zeta$ and $\phi_h(\zeta) = 1 + \beta_h \zeta$, the gradient Richardson number is:  
  
$$Ri_g = \frac{\zeta \phi_h(\zeta)}{\phi_m(\zeta)^2} = \frac{\zeta(1 + \beta_h \zeta)}{(1 + \beta_m \zeta)^2}$$  
Rearranging into standard quadratic form $A \zeta^2 + B \zeta + C = 0$:  
  
$$(\beta_m^2 Ri_g - \beta_h) \zeta^2 + (2 \beta_m Ri_g - 1) \zeta + Ri_g = 0$$  
The discriminant simplifies analytically as the $4\beta_m^2 Ri_g^2$ terms cancel out:  
  
$$\Delta = (2 \beta_m Ri_g - 1)^2 - 4(\beta_m^2 Ri_g - \beta_h) Ri_g = 1 + 4(\beta_h - \beta_m) Ri_g$$  
The physically meaningful root ($\zeta \ge 0$ for $Ri_g \ge 0$) is:  
  
$$\zeta(Ri_g) = \frac{1 - 2\beta_m Ri_g - \sqrt{1 + 4(\beta_h - \beta_m) Ri_g}}{2(\beta_m^2 Ri_g - \beta_h)}$$  
Once $\zeta(Ri_g)$ is evaluated, the closures and emergent Prandtl number are computed directly without approximations:  
  
$$f_m(Ri_g) = \frac{1}{\phi_m(\zeta)^2}, \quad f_h(Ri_g) = \frac{1}{\phi_m(\zeta) \phi_h(\zeta)} \implies Pr_t(Ri_g) = \frac{f_m}{f_h} = \frac{\phi_h(\zeta)}{\phi_m(\zeta)}$$  
## 2. Refactored Module Architecture (src/Closures.jl)  
Julia  
  
module Closures  
  
using LinearAlgebra  
  
export AbstractStabilityClosure,  
       MOSTEqualSlope,  
       MOSTGeneral,  
       GrachevVSBL,  
       GeometricWrapper,  
       ClosureState,  
       evaluate_closure,  
       compute_diffusivities!  
  
const DEFAULT_SHEAR2_MIN = 1e-12  
  
# -------------------------------------------------------------------  
# Stage 1: Closure State & Type Abstractions  
# -------------------------------------------------------------------  
  
"""  
    ClosureState{T}  
  
Encapsulates the non-dimensional outputs of a stability closure evaluation.  
- `fm`: Momentum stability function  
- `fh`: Heat stability function  
- `Prt`: Emergent turbulent Prandtl number (fm / fh)  
"""  
struct ClosureState{T<:Real}  
    fm::T  
    fh::T  
    Prt::T  
end  
  
abstract type AbstractStabilityClosure end  
  
"""  
    MOSTEqualSlope(beta=5.0)  
  
Exact analytical closure for equal-slope linear MOST (β_m = β_h = β).  
Inverts directly to f_m = f_h = (1 - β * Ri_g)^2 with critical Richardson number Ri_c = 1 / β.  
"""  
struct MOSTEqualSlope <: AbstractStabilityClosure  
    beta::Float64  
    Ri_c::Float64  
  
    MOSTEqualSlope(beta=5.0) = new(beta, 1.0 / beta)  
end  
  
"""  
    MOSTGeneral(beta_m=5.0, beta_h=5.0)  
  
General linear MOST closure for unequal slopes (β_m ≠ β_h).  
Solves the exact quadratic relation ζ(Ri_g) to compute φ_m(ζ) and φ_h(ζ).  
"""  
struct MOSTGeneral <: AbstractStabilityClosure  
    beta_m::Float64  
    beta_h::Float64  
    Ri_c::Float64  
  
    MOSTGeneral(beta_m=5.0, beta_h=5.0) = new(beta_m, beta_h, 1.0 / beta_m)  
end  
  
"""  
    GrachevVSBL()  
  
SHEBA / Grachev et al. (2007) non-linear stability functions for Very Stable Boundary Layers (VSBL).  
"""  
struct GrachevVSBL <: AbstractStabilityClosure end  
  
# -------------------------------------------------------------------  
# Stage 2: Closure Evaluators returning ClosureState  
# -------------------------------------------------------------------  
  
@inline function evaluate_closure(c::MOSTEqualSlope, Ri_g::Float64)::ClosureState{Float64}  
    if Ri_g >= c.Ri_c  
        return ClosureState(0.0, 0.0, 1.0)  
    end  
    f = (1.0 - c.beta * Ri_g)^2  
    return ClosureState(f, f, 1.0)  
end  
  
@inline function evaluate_closure(c::MOSTGeneral, Ri_g::Float64)::ClosureState{Float64}  
    if Ri_g >= c.Ri_c  
        return ClosureState(0.0, 0.0, 1.0)  
    end  
  
    # Solve quadratic equation for ζ(Ri_g)  
    bm = c.beta_m  
    bh = c.beta_h  
      
    A = bm^2 * Ri_g - bh  
    B = 2.0 * bm * Ri_g - 1.0  
    discriminant = max(0.0, 1.0 + 4.0 * (bh - bm) * Ri_g)  
      
    # Evaluate root  
    zeta = abs(A) > 1e-12 ? (1.0 - 2.0 * bm * Ri_g - sqrt(discriminant)) / (2.0 * A) : Ri_g  
  
    phi_m = 1.0 + bm * zeta  
    phi_h = 1.0 + bh * zeta  
  
    fm = 1.0 / (phi_m^2)  
    fh = 1.0 / (phi_m * phi_h)  
    Prt = phi_h / phi_m  
  
    return ClosureState(fm, fh, Prt)  
end  
  
@inline function evaluate_closure(::GrachevVSBL, Ri_g::Float64)::ClosureState{Float64}  
    fm = (1.0 + 5.0 * Ri_g)^(-0.8)  
    fh = (1.0 + 5.0 * Ri_g)^(-2.0)  
    Prt = fm / fh  
    return ClosureState(fm, fh, Prt)  
end  
  
# -------------------------------------------------------------------  
# Stage 3: Topological Geometric Correction Layer  
# -------------------------------------------------------------------  
  
struct GeometricWrapper  
    delta_0::Float64  
    enabled::Bool  
end  
  
GeometricWrapper(delta_0=1e-3) = GeometricWrapper(delta_0, true)  
  
@inline function evaluate_manifold_correction(gw::GeometricWrapper, d_fold::Float64)::Float64  
    return gw.enabled ? tanh(d_fold / gw.delta_0) : 1.0  
end  
  
# -------------------------------------------------------------------  
# Stage 4 & 5: Column Execution Pipeline  
# -------------------------------------------------------------------  
  
"""  
    compute_diffusivities!(Km, Kh, Pr_t, z, u, v, theta, closure, geom_wrapper, d_fold_profile;   
                           shear2_min=DEFAULT_SHEAR2_MIN, ka=0.4, l0=30.0)  
  
Executes the 3-tiered closure pipeline across an atmospheric column.  
"""  
function compute_diffusivities!(  
    Km::Vector{Float64},  
    Kh::Vector{Float64},  
    Pr_t::Vector{Float64},  
    z::Vector{Float64},  
    u::Vector{Float64},  
    v::Vector{Float64},  
    theta::Vector{Float64},  
    closure::AbstractStabilityClosure,  
    geom_wrapper::GeometricWrapper,  
    d_fold_profile::Vector{Float64};  
    shear2_min::Float64=DEFAULT_SHEAR2_MIN,  
    ka::Float64=0.4,  
    l0::Float64=30.0,  
    g::Float64=9.81  
)  
    Nz = length(z)  
    theta0 = theta[1]  
  
    for i in 2:(Nz - 1)  
        # 1. Gradients with explicit regularization  
        dz = z[i+1] - z[i-1]  
        du = (u[i+1] - u[i-1]) / dz  
        dv = (v[i+1] - v[i-1]) / dz  
        dth = (theta[i+1] - theta[i-1]) / dz  
  
        shear_sq = max(du^2 + dv^2, shear2_min)  
        shear = sqrt(shear_sq)  
  
        # 2. Mixing Length (Blackadar)  
        l_mix = (ka * z[i]) / (1.0 + (ka * z[i]) / l0)  
  
        # 3. Richardson Number  
        N2 = (g / theta0) * dth  
        Ri_g = max(0.0, N2 / shear_sq)  
  
        # 4. Local Physics: ClosureState Evaluation  
        state = evaluate_closure(closure, Ri_g)  
  
        # 5. Base Diffusivities  
        base_K = (l_mix^2) * shear  
        Km_base = base_K * state.fm  
        Kh_base = base_K * state.fh  
  
        # 6. Global Geometry: Topological Suppression Factor C(d_fold)  
        C_fold = evaluate_manifold_correction(geom_wrapper, d_fold_profile[i])  
  
        # Final Column Quantities  
        Km[i] = Km_base * C_fold  
        Kh[i] = Kh_base * C_fold  
          
        # Pr_t remains invariant under amplitude scaling C(d_fold)  
        Pr_t[i] = state.Prt  
    end  
  
    # Boundary conditions  
    Km[1], Kh[1], Pr_t[1] = Km[2], Kh[2], Pr_t[2]  
    Km[Nz], Kh[Nz], Pr_t[Nz] = Km[Nz-1], Kh[Nz-1], Pr_t[Nz-1]  
  
    return nothing  
end  
  
end # module Closures  
## Benefits of the Refined Design  
1. **Exact Mathematical Separation:** MOSTEqualSlope handles the direct polynomial $f(Ri_g) = (1 - \beta Ri_g)^2$, while MOSTGeneral solves the full quadratic $\zeta(Ri_g)$ relation when $\beta_m \neq \beta_h$.   
2. **ClosureState{T} Encapsulation:** Returning a unified struct avoids returning loose tuples and ensures $Pr_t = f_m / f_h$ is recorded directly alongside $f_m$ and $f_h$.   
3. **Invariance of Transport Anisotropy:** Because $C(d_{\text{fold}})$ multiplies both $K_m$ and $K_h$ equally, $Pr_t[i] = K_m / K_h$ remains mathematically identical to state.Prt, ensuring slow-manifold geometry regulates turbulence amplitude without corrupting flux ratios.   
4. **Configurable Regularization:** Floor parameters like shear2_min are explicitly named and exposed rather than hardcoded inline.  
