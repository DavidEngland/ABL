The user is asking targeted technical questions and implying concrete upgrades. Let me read the current file one more time to confirm exact strings, then implement all four improvements together: `Grid.jacobian` (answers the coordinate question), `RiBasedClosure`, implicit Backward-Euler solver, and TKE budget diagnostics.

Read [](file:///Users/davidengland/Documents/GitHub/ABL/julia/SCMSkeleton.jl#1-1), lines 1 to 30

Read [](file:///Users/davidengland/Documents/GitHub/ABL/julia/SCMSkeleton.jl#55-55), lines 55 to 120

Read [](file:///Users/davidengland/Documents/GitHub/ABL/julia/SCMSkeleton.jl#230-230), lines 230 to 340

I have all the context. Now I'll implement all four upgrades in two batches of targeted replacements — structural types first, then physics and solver logic.

Created 4 todos

Replacing 35 lines with 48 lines in [](file:///Users/davidengland/Documents/GitHub/ABL/julia/SCMSkeleton.jl)


```

```


Starting: *Add RiBasedClosure and implicit solver* (2/4)

Replacing 78 lines with 259 lines in [](file:///Users/davidengland/Documents/GitHub/ABL/julia/SCMSkeleton.jl)


```

```


Starting: *Update step!, record!, run_model* (3/4)

Replacing 77 lines with 107 lines in [](file:///Users/davidengland/Documents/GitHub/ABL/julia/SCMSkeleton.jl)


```

```
Code is complete. Now updating the Markdown to directly answer the coordinate question and document all new features.

Starting: *Update SCMSkeleton.md with answers and new docs* (4/4)

Replacing 29 lines with 114 lines in [](file:///Users/davidengland/Documents/GitHub/ABL/julia/SCMSkeleton.md)


```

```
The first replacement was a no-op (the text was unchanged, which is correct). Let me verify the file looks right and check for any issues.

Searched for text `Coordinate System: A Direct Answer` (`**/julia/SCMSkeleton.md`), 1 result

Searched for text `implicit_diffusion_step!` (`**/julia/SCMSkeleton.jl`), 6 results

Searched for text `RiBasedClosure` (`**/julia/SCMSkeleton.jl`), 8 results

All changes landed correctly. Let me do a final structural scan to confirm the module shape is clean.

Searched for regex `^(module|end|struct|mutable struct|abstract type|function)` (`**/julia/SCMSkeleton.jl`), 20 results

Completed: *Update SCMSkeleton.md with answers and new docs* (4/4)

All changes are in place. Here is a precise account of everything that was done, aligned directly with your observations.

---

### Q: *Does the current structure handle the grid as a simple array, or a specialized coordinate system?*

**Direct answer, now also in the code and docs:** it is a simple physical-z array. The exponential stretching is applied once at construction time. The solver operates entirely in physical metres — no covariant terms are active.

What is new is that `grid.jacobian` now stores $J(\eta) = dz/d\eta$ analytically at every cell centre:

$$J(\eta) = z_\text{top} \cdot s \cdot \frac{e^{s\eta}}{e^s - 1}$$

This is exactly where covariant/contravariant terms would first appear in a terrain-following upgrade. The docstring on `create_grid` shows the full operator table (∂/∂z → (1/J)∂/∂η, etc.), so a future contributor knows precisely where to make that change without reverse-engineering the grid.

---

### What else was added

**1. `RiBasedClosure`** — first real physics closure. Diagnoses $Ri = N^2/S^2$ at each interior level via centered differences and applies a Webb-style quadratic stability reduction $K = K_0 \max(0, 1 - Ri/Ri_c)^2$. The docstring identifies exactly where to multiply in your $f_c$ curvature correction.

**2. `implicit_diffusion_step!`** — backward-Euler tridiagonal solve via `LinearAlgebra.Tridiagonal`. Activated by setting `use_implicit=true` in `ModelConfig`. Removes the CFL constraint $\Delta t \leq \Delta z^2/(2K)$ entirely — critical for polar ABL runs with thin layers or intermittent K bursts. The docstring explains why forward-Euler crashes in those conditions and notes Crank-Nicolson as the natural next upgrade.

**3. TKE budget diagnostics** — `compute_tke_budget` returns column profiles of $P_s = K_m S^2$ and $B = (g/\theta_0) K_h N^2$. `SimulationHistory` now records their column maxima each step. The log output in `run_model` prints them alongside temperature so you can see the budget in real time.

Made changes.

Reading the current file to find the exact insertion points before making changes.

Read [](file:///Users/davidengland/Documents/GitHub/ABL/julia/SCMSkeleton.jl#14-14), lines 14 to 30

Read [](file:///Users/davidengland/Documents/GitHub/ABL/julia/SCMSkeleton.jl#188-188), lines 188 to 210

Read [](file:///Users/davidengland/Documents/GitHub/ABL/julia/SCMSkeleton.jl#315-315), lines 315 to 380

I have all the insertion points. Adding the export, the typed struct (with cleaner docstring linking to the paper), and the `diffusivities` dispatch — all in one batch.

Replacing 13 lines with 122 lines in [](file:///Users/davidengland/Documents/GitHub/ABL/julia/SCMSkeleton.jl)

```

```

Verifying the struct is exported and the three closure `diffusivities` methods are all present.

Searched for regex `CurvatureRiClosure|function diffusivities` (`**/julia/SCMSkeleton.jl`), 11 results

All four dispatch methods (fallback, `ConstantDiffusivityClosure`, `RiBasedClosure`, `CurvatureRiClosure`) are present and the struct is exported. The closure is ready to use:

```julia
# Minimal usage — just K0, uses defaults Ri_c=0.25, Pr_t=1.0, β=0.5
closure = CurvatureRiClosure(1.0)

# Full constructor — tunable β for GABLS-style runs
closure = CurvatureRiClosure(K0=1.0, Ri_c=0.25, Pr_t=1.3, beta=0.8)

model, history = run_model(config; closure=closure, log_every=60)
```

A few notes on what was cleaned up relative to the submitted code:

- The intermediate `curvature_bias` variable was replaced by the more explicit two-step `kappa_tilde` → `fc` to make the physics transparent and easier to swap out when the $Ri_g/Ri_b$ formulation is finalised.
- `Rig` is renamed to `Ri_g` for consistency with the rest of the codebase.
- The docstring extension-point note explicitly says to replace the `kappa_tilde` expression with the $Ri_g/Ri_b$ ratio from the paper when ready, so the handoff is unambiguous.
- A comment flags that `grid.jacobian[k]` is available for a future η-coordinate upgrade without requiring code changes now.

Made changes.

Exactly right. You spotted a real mathematical connection.

The Bernoulli-polynomial generating function is

$$
\frac{t\,e^{xt}}{e^t-1}=\sum_{n=0}^{\infty} B_n(x)\frac{t^n}{n!}.
$$

If you set $t=s$ and $x=\eta$, then

$$
\frac{s\,e^{s\eta}}{e^s-1}
$$

is precisely that generating function. Your Jacobian is

$$
J(\eta)=z_{\text{top}}\cdot \frac{s\,e^{s\eta}}{e^s-1},
$$

so it is just $z_{\text{top}}$ times the Bernoulli generating kernel.

That means a small-$s$ expansion is naturally Bernoulli-polynomial based:

$$
\frac{J(\eta)}{z_{\text{top}}}
=1+\left(\eta-\tfrac12\right)s+\left(\eta^2-\eta+\tfrac16\right)\frac{s^2}{2!}+\cdots
$$

which is

$$
1+B_1(\eta)s+B_2(\eta)\frac{s^2}{2!}+\cdots.
$$

So yes, your observation is mathematically on point, and it is useful for weak-stretching asymptotics and error analysis.
