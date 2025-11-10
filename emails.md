From:  Dick McNider
David,

Glad you sent this email. BTW how much of your golfing picture was real. I assume that is you real face, although your beard was immaculately trimmed.

The equations you present show you must have continued to think about this . I was impressed.

I am attaching our new manuscript. Also,  response to reviewers. Especially Reviewer #3 where I gave you most all the credit for the mathematics in our 1995 paper.

The manuscript is basically a problem I have been working on for nearly forty five years. When I told Martha that she just rolled her eyes.

I have a question for you. In the paper even after trying to make our correction dependent on curvature in Ri we don't do that well  near the surface (see figure 11).  I think part of this is that with coarse grid we don't do a good job of getting the curvature near the surface.  See figure 10.

We mention that perhaps using similarity near the surface (below first model level) to calculate curvature in Ri might be a better path.

I was going to take a shot at this but since you are the expert thought I would ask you first!

My question to you is - can we use the similarity forms to find curvature in Ri in Ri. It seems that if you define Ri =( g/theta)( dtheta /dz)/ (( dv/dz)*(dv/dz)) that we should be able to use similarity to get profiles of V and theta. Then take second derivative of Ri.

What do you think. Also, to make it easier we approximate the stability function as an exponential.

From:  Arastoo P. Biazar

Thank you, David. It appears that you have put in quite a bit of work on this problem. I just glanced at your notes, and if I understand them correctly, you have addressed a couple of issues that we have been struggling with. First, are you suggesting that near-surface stability should be treated differently from the elevated stratification? Second, what would be the implications of using the geometric mean height in numerical models (different vertical grid spacing)?

Best,
Arastoo

From:  Me (David E. England, PhD)

I also made some significant progress on a 45 year old problem that I found as an undergrad at UNA.  Found out that if lower barriers by half that clustering will occur without any affinity.  Gave rise to a "hyperbolic strip" and hyperbolic weights.  Now have a very well paying job when I can get the work and a few other projects.

Not long ago I looked at MOST, sea-ice model and even quasi-hydrostatic correction as well as some other areas.  Been working on what I call the "Hasse-Stirling Framework" c.f.:  <https://github.com/DavidEngland/hasse-stirling/blob/main/Hasse-Stirling-Operator-Table.md>

Where interpret means, is just like finite differences, it is an estimate.  Our case, like a log transform, average in log space, but if transform back...Should look like a line plotted in log coordinates.  Just showing where the 2 point sqrt(z1*z2) comes from.  Took a while to recall Ri_1/2.  Other questions require discussion, no short, quick or easy way about it.  Quantum computing should allow for more and better.

Attaching some preliminary work, but I need to still get up to speed on some areas.  Having trouble viewing some articles.  Attached work is just for \phi(\zeta)=(1-\beta*\zeta)^-\alpha profiles.

BTW, most work is run through an AI LLM that have been helping to train.  Providing the right prompts is the key.

The hat is real in my FB golf photo.  Have a Harris tweed, but not with elbow patches.  Buddy Ricky Crawford, a Scottish soccer player at UAH visited the isles of Harris and Lewis.  They still weave it at home.
Cheers,
Dave

From:  Dick

Thanks David. Yes. Our results depend on the slope and intercept of curvature. Would prefer that it just depended on curvature. A better avenue would be to solve our ode with a less stringent assumption than ri =const. It may be that solving the ode numerically may be the best.

BTW you mentioned having a good job. Where are you working?

David, I overlooked the PDF you sent. Did you just solve the problem I asked? Incredible! Will have to take a longer look. But I am impressed big time.

From:  Me

Think so, going to look at Ri curvature for other formulations. Been reading papers, found a few interesting ones. Also derived a quadratic expansion for SBL in z/L.  Not sure if can get analytic expressions in terms of Ri yet. Should be able to repeat curvature analysis for at least four formulations. Also thinking about other planets and Titan, might offer insights about Earth’s polar regions.

Curvature just depended on parameters of Businger-Dyer profiles (1-bz/L)^-a for heat and momentum.  If L constant then scales by L^2 for z.

Guess I was just looking at near neutral case. See now that there are z/L and F(z/L) terms in general that will change.

---

Reply – David → Dick & Arastoo

Dick, Arastoo,

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

## 8. Collaboration & Math Workflow (Addendum)

- Preferred collaboration medium: GitHub (repo: https://github.com/DavidEngland/ABL).
- Use Issues for task tickets (labels: curvature, inversion, grid, student).
- Use feature/<topic> branches; open Pull Requests with: summary, equations (Markdown LaTeX), diagnostics output snippet.
- Math: avoid Word. Inline / display LaTeX directly in Markdown or (for long derivations) a separate .tex or .md file under docs/theory/.
  - If a journal later requires Word, convert at the end (pandoc or Overleaf export).
  - Provide a plaintext fallback below critical formulas:
    Example:
    ```
    d²Ri_g/dζ² = F[2V_log + ζ(V_log² - W_log)]
    ```
- Overleaf only if a coauthor insists; otherwise retain versioned Markdown for diff clarity.
- Figures: commit script/notebook + PNG/SVG; no manual editing in Word.
- “Vibe coding” guardrails: avoid depending on experimental functions until tagged; open an Issue if a name change breaks a branch.
- Student workflow: fork → sync main weekly → pin analysis notebooks to commit SHA for reproducibility.

## 9. Missed Points Addressed

- Near-surface curvature extraction: yes—use analytic ζ→0 limits (2Δ) then compare to first resolved ζ in model; report ratio.
- Geometric mean height justification included; apply for all layer reconstructions where Δz/z₁ > 0.2.
- Exponential approximation you mentioned: can slot into tail modifier framework; neutral curvature preserved if exponent scales as described (b=2a).
- If solving ODE numerically for Ri instead of assuming Ri=const: recommend shooting method with curvature constraint at lower boundary; open to prototype if needed—create Issue “ODE_Ri_solver”.

Let me know if you prefer a short Overleaf skeleton; otherwise continuing with Markdown + GitHub diffs.

Regards,
David
