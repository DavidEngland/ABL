"""
    CurvatureRiClosure(K0, Ri_c, beta)

An extension of the Ri-based closure that applies a grid-dependent
correction factor (fc) based on the curvature of the temperature profile.

Fields
──────
- `K0`:   Neutral diffusivity (m²/s).
- `Ri_c`: Critical Richardson number.
- `beta`: Sensitivity parameter for the curvature correction.
"""
struct CurvatureRiClosure <: AbstractClosure
    K0::Float64
    Ri_c::Float64
    Pr_t::Float64
    beta::Float64  # Scaling factor for curvature correction
end

# Default constructor
CurvatureRiClosure(K0::Float64) = CurvatureRiClosure(K0, 0.25, 1.0, 0.5)

function diffusivities(closure::CurvatureRiClosure, model::SCMModel)
    grid  = model.grid
    state = model.state
    nz    = grid.nz
    g     = 9.81
    theta_ref = sum(state.theta) / nz

    km = fill(1.0e-4, nz)
    kh = fill(1.0e-4, nz)

    for k in 2:nz-1
        # 1. Standard First Derivatives (Central Difference)
        dz_span = grid.z[k+1] - grid.z[k-1]
        d_theta = state.theta[k+1] - state.theta[k-1]
        du      = state.u[k+1]     - state.u[k-1]
        dv      = state.v[k+1]     - state.v[k-1]

        N2 = (g / theta_ref) * d_theta / dz_span
        S2 = (du^2 + dv^2) / dz_span^2
        Rig = S2 > 1.0e-9 ? N2 / S2 : 1.0e6

        # 2. Curvature Diagnosis (Second Derivative of Theta)
        # Using the Jacobian (J) here ensures the curvature is
        # consistent with the stretched grid geometry.
        dz_up = grid.z[k+1] - grid.z[k]
        dz_dn = grid.z[k]   - grid.z[k-1]
        d2_theta = ((state.theta[k+1] - state.theta[k])/dz_up -
                    (state.theta[k] - state.theta[k-1])/dz_dn) / grid.dz[k]

        # 3. Compute Correction Factor fc
        # If curvature is positive (stable inversion sharpening),
        # fc enhances mixing to prevent numerical decoupling.
        # This is a placeholder for your specific Ri_g/Ri_b logic.
        curvature_bias = 1.0 + closure.beta * (d2_theta * grid.z[k] / max(d_theta/dz_span, 1e-4))
        fc = max(0.1, min(2.0, curvature_bias))

        # 4. Apply Stability Function + Correction
        if Rig < closure.Ri_c
            f_stability = (1.0 - Rig / closure.Ri_c)^2
            km[k] = closure.K0 * f_stability * fc
            kh[k] = closure.K0 * f_stability * fc / closure.Pr_t
        end
    end
    return km, kh
end

Here’s the good news, David: what you’ve written is already very close to a publishable prototype of a curvature‑aware Ri closure. But there are a few places where the physics can be sharpened, the numerics made safer, and the conceptual framing aligned with what you’ve been building—namely:

• a dynamic, state‑dependent critical Ri,
• derived from higher‑order turbulence structure,
• with curvature acting as a proxy for local intermittency and inversion sharpening,
• while still remaining compatible with existing Ri‑based schemes McNider wants to preserve.


Below is a structured critique and refinement path that keeps your intent intact but strengthens the formulation.

---

🔍 What You’ve Built (in conceptual terms)

Your closure currently does three things:

1. Diagnoses local Ri using first derivatives

This is standard and fine, though the central difference over stretched grids can be noisy.

2. Diagnoses curvature using a second derivative of θ

This is the interesting part. The curvature term is acting as a proxy for inversion sharpening, which is exactly where MOST and classical Ri closures fail.

3. Modulates mixing via a multiplicative factor `fc`

This is your “dynamic Ri_c” in disguise:
you’re effectively saying

Ri_{c,\mathrm{eff}} = \frac{Ri_c}{f_c}


which is a perfectly reasonable way to embed higher‑order structure into a first‑order closure.

---

⚙️ Where the Physics Can Be Strengthened

1. The curvature ratio is dimensionally fragile

Right now you compute:

\text{curvature\_bias} = 1 + \beta \frac{d^2\theta}{dz^2} \frac{z}{d\theta/dz}


This ratio is physically interpretable (curvature relative to slope), but numerically it will blow up whenever the gradient is small.

You already guard with max(d_theta/dz_span, 1e-4), but that’s a blunt instrument.

Better approach: use a nondimensional curvature parameter

Define:

C = \frac{z}{\theta} \frac{d^2\theta}{dz^2}


or even better:

C = \frac{d^2\theta/dz^2}{|d\theta/dz|/L}


where \(L\) is a local mixing length (grid‑scale or Blackadar‑type).

This keeps the curvature measure tied to the local turbulence scale, not the raw grid spacing.

---

2. The correction factor should be asymmetric

Right now positive and negative curvature are treated symmetrically except for the sign.

But physically:

• Positive curvature (inversion sharpening) → enhance mixing
• Negative curvature (convective curvature) → do not suppress mixing too much


A better formulation is:

f_c = 1 + \beta \max(C, 0)


and then clamp to \([0.1, 2.0]\).

This prevents pathological suppression of mixing in weakly unstable layers.

---

3. The stability function should incorporate fc inside the Ri argument

Right now you do:

K \propto f_{\text{stability}} \cdot f_c


But the physics says:

• curvature modifies the effective critical Ri
• not the diffusivity directly


So instead use:

Ri_{c,\mathrm{eff}} = Ri_c \, f_c


and then:

f_{\text{stability}} = \left(1 - \frac{Ri}{Ri_{c,\mathrm{eff}}}\right)^2


This is much more interpretable and aligns with your “dynamic critical Ri” idea.

---

4. Consider using the buoyancy flux sign to gate curvature effects

Curvature only matters when the layer is stably stratified.

So add:

if N2 > 0
    # apply curvature correction
else
    fc = 1.0
end


This prevents curvature from interfering with unstable layers.

---

🧠 A More Physically Grounded Version of Your Correction

Here’s a refined conceptual formula you might like:

C = \frac{d^2\theta/dz^2}{|d\theta/dz|/L}


f_c = 1 + \beta \max(C, 0)


Ri_{c,\mathrm{eff}} = Ri_c \, f_c


K = K_0 \left(1 - \frac{Ri}{Ri_{c,\mathrm{eff}}}\right)^2


This is clean, interpretable, and numerically stable.

---

🧩 What You Have Is Already a Prototype of a Higher‑Order Closure

Your curvature term is essentially a second‑order structure function embedded inside a first‑order closure. That’s exactly the right direction for:

• dynamic Ri_c
• intermittency-aware mixing
• polar inversion physics
• ML emulation later


If you want, I can help you:

• rewrite the closure in a more compact, physically transparent form
• derive the nondimensional curvature parameter cleanly
• embed this into a dynamic Ri_c formulation
• compare it to Mellor–Yamada level‑2.5 and QNSE
• or even sketch how ML would emulate the curvature term without breaking physics


Just tell me which direction you want to push next.
