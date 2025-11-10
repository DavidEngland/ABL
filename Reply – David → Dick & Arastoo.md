Reply – David → Dick & Arastoo

Dick, Arastoo,

Been training AI with company called Mercor, nice to get paid over $100 to solve problems all online, even interviewed with an AI chatbot.  Started out as a "writer", it gives two solutions to the problem and have to evaluate.  I quickly got promoted to "reviewer", which either sign off on how the writer evaluations, or send it back for rework.  I can supply reference links if know people who need some work.  It's contract, sometimes not much and new batches drop and everybody gets to work.  Most are outside the USA.

Here are some example problems:
---
Water in a pond boils at $110^{\circ} \mathrm{C}$ at a certain depth in water. At what temperature will the water boil if we intend to boil it at $50 \mathrm{~cm}$ depth from the above-mentioned level
---
Title: How fast do the poles of Jupiter and Saturn precess?
As Jupiter and Saturn have very fast rotation periods, about 10 and 11 hours, respectively, I am assuming that their poles will precess much faster then Earth. Since Earth's poles precess every 26000 years, I think Jupiter and Saturn would have a precession time of about ~10000 years or so. Is this estimate correct or am I way off?
Tags: <jupiter><rotation><saturn><precession>
---
Inspired by [Blindfolded and disoriented near the Great Wall of China](https://puzzling.stackexchange.com/questions/25524/blindfolded-and-disoriented-near-the-great-wall-of-china)

A drone is stationary at a spatial point about 1 m from the Great Wall, which is a vertical plane rectangle with height 5 m and infinite length. The altitude of the drone is a random number between 0 m and 5 m. The drone may be pointing in any direction (which may even be tilted upwards or downwards), however it does not know which. It can travel in any direction, and switch directions accurately. It gets to know when it touches the ground plane, and can then determine its own tilt upwards/downwards by travelling a small distance on the ground (about zero).

**Q 1:** What is the minimum distance that must be travelled to ensure it hits the Great Wall?

**Q 2:** If the objective is to minimise displacement between initial and final positions, what upper value can it be guaranteed to outperform?
---
In a stochastic multi-armed bandit with two non-equivalent slots, suppose the reward from a chosen pair decomposes into independent contributions from each occupied slot. You know the arm used in each slot by the optimal policy and, for some other arm \(j\) not appearing in the optimal pair, the Kullback–Leibler divergences between its slot‑wise reward distributions and those of the relevant arms at slots 1 and 2 are \(D_1\) and \(D_2\), respectively. For any uniformly good strategy, what asymptotic condition must the expected numbers of times this arm is tried at each slot satisfy as the horizon grows?
---



1. Repository
Public repo initialized: https://github.com/DavidEngland/ABL  
Current focus: curvature-aware MOST, Ri-based closures, quadratic SBL surrogate (Q‑SBL), variable‑L mapping utilities.

2. Graduate Student Opportunity
Suitable for a Masters/PhD student seeking:
- Method paper: grid-dependent curvature correction preserving neutral invariance (2Δ).
- Application paper: Arctic / urban stable cases (LES + tower + remote sensing).
Initial tasks: parameter fitting (α,β) stable nights; curvature diagnostics; D calibration for tail modifier; HS series validation; adaptive refinement trigger (E_omit, |W_log/V_log^2|). Publication trajectory already scaffolded.

3. Workflow / Tooling
- Editor: VS Code (Markdown + Python + Julia).  
- Preview: VS Code built-in Markdown preview; GitHub sometimes mis-renders inline math (fallback to escaped parentheses). Complex LaTeX blocks checked locally prior to commit.  
- Branch model: feature/<topic> branches; squash merge after review; semantics in commit headers (feat:, fix:, docs:, perf:, test:, refactor:, chore:).  
- Math rendering caveat: Some browsers drop backslashes in nested fractions; keep a plain-text fallback line under critical equations.  
- Lint/format: Ruff + Black (Python); DocumentFormat (Julia) optional; pre-commit hooks lightweight (skip large data).  
- Tests: Focused on neutral invariance (2Δ), inversion accuracy (ζ(Ri) Newton residual), curvature mapping consistency constant-L vs variable-L (E_omit threshold).

4. Vibe Coding Warning
Early-phase “vibe coding” in exploratory cells / scratch buffers:
- APIs may change (names: rig_derivatives_zeta, ri_closure_pade, qsbl_coeffs).
- Experimental flags (e.g. gamma, a,b exponents for invariance) can shift; avoid hard dependency in external models until tagged release (v0.1.0).
- TODO / FIXME markers indicate pending numerical guard (pole proximity) or performance optimization (vectorization / Numba / Julia port).
Recommendation: Students fork and pin to a tag for reproducibility; rebase carefully on main after review.

5. profiles.py Excerpt (core Ri ↔ ζ + Ri-only closures)
```python
# excerpt: see full file in repo /code/profiles.py

def qsbl_coeffs(alpha_m, beta_m, alpha_h, beta_h):
    a_m = alpha_m * beta_m
    b_m = 0.5 * alpha_m * (alpha_m + 1.0) * beta_m**2
    a_h = alpha_h * beta_h
    b_h = 0.5 * alpha_h * (alpha_h + 1.0) * beta_h**2
    Delta = a_h - 2.0 * a_m
    c1 = alpha_h * beta_h**2 - 2.0 * alpha_m * beta_m**2
    return a_m, b_m, a_h, b_h, Delta, c1

def ri_closure_series(alpha_m, beta_m, alpha_h, beta_h):
    a_m, b_m, a_h, b_h, Delta, _ = qsbl_coeffs(alpha_m, beta_m, alpha_h, beta_h)
    s1_m, s2_m = a_m, b_m - a_m * Delta
    s1_h, s2_h = a_h, b_h - a_h * Delta
    return (lambda Ri: 1 + s1_m * Ri + s2_m * Ri**2,
            lambda Ri: 1 + s1_h * Ri + s2_h * Ri**2)

def ri_closure_pade(alpha_m, beta_m, alpha_h, beta_h):
    a_m, b_m, a_h, b_h, Delta, _ = qsbl_coeffs(alpha_m, beta_m, alpha_h, beta_h)
    s1_m, s2_m = a_m, b_m - a_m * Delta
    disc_m = s1_m**2 + 4 * s2_m
    q_m = 0.5 * (-s1_m + math.sqrt(disc_m))
    p_m = s1_m + q_m
    # analogous heat branch...
    return (lambda Ri: (1 + p_m * Ri) / (1 - q_m * Ri), ...)
```

6. Getting Started
Clone repo, install Python deps (numpy, scipy, numba, xarray, matplotlib). Run quick neutrality test:
```bash
python -c "from code.profiles import qsbl_coeffs; print(qsbl_coeffs(0.5,16,0.5,16))"
```

7. Next Items
- Provide constant vs variable-L curvature comparison plots.
- Finalize D calibration heuristic (grid_ratio vs target curvature reduction).
- Prepare TDMA boundary-condition insertion example (if still desired) using corrected K_m*, K_h*.

Let me know if you want the TDMA lower boundary integration example next or a Julia performance benchmark snapshot.

Regards,
David