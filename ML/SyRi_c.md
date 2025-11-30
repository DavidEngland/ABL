# Symbolic regression workflow for dynamic critical Richardson number

A compact, physics-informed expression for \(\mathrm{Ri_c^*}\) can sharpen turbulence onset/cessation timing while remaining safe and interpretable. Below is a concrete, end-to-end recipe to implement, validate, and deploy symbolic regression for \(\mathrm{Ri_c^*} = f(\text{state})\) with operational guard rails.

---

## Objective and success metrics

- **Primary objective:**
  \(\mathrm{Ri_c^*}\) predicts collapse/onset timing from state, improving event skill versus fixed \(\mathrm{Ri_c}\) (e.g., 0.25).

- **Output bounds and behavior:**
  \(\mathrm{Ri_c^*} \in [0.15, 2.0]\); monotone with respect to select drivers where physics demands (e.g., increase with stratification, decrease with shear).

- **Regression metrics:**
  **MAE** and **RMSE** on \(\mathrm{Ri_c^*}\) at held sites/times.

- **Event metrics:**
  **Precision**, **recall**, and **F1** for predicted collapse/restore within \(\pm \Delta t\) (e.g., 10 min).

- **Operational criteria:**
  **Closed-loop stability** (no solver crashes), **energy consistency** (no spurious sources/sinks), and **surface bias change** within tolerance.

- **Example success thresholds:**
  **MAE** < 0.05 and **F1** increase ≥ 20% vs baseline fixed \(\mathrm{Ri_c}\).

---

## Data preparation and splits

- **Sources:**
  LES ensembles (e.g., GABLS-style), tower observations (SHEBA, CASES-99, ARM), synthetic MOST sweeps.

- **Labels and targets:**
  **Turbulence ON/OFF** via TKE/flux thresholds; derive per-event \(\mathrm{Ri_c^*}\) as the \(\mathrm{Ri}\) where state flips in a controlled reanalysis.

- **Sampling and folds:**
  **Site/time-based splits** to avoid leakage; stratify by regime (stable, neutral, weak shear) and by \(\Delta z\).

- **Quality controls:**
  **Outlier filters** (sensor spikes, unphysical gradients), **consistency checks** (signs of \(\Gamma = d\theta/dz\), \(S = |\partial U/\partial z|\)), and **minimum data density** per window.

---

## Feature engineering and scalings

- **Physics-first features:**
  \(\Gamma = d\theta/dz\), \(S = |\partial U/\partial z|\), \(\zeta = z/L\), \(\Delta z\), \(z_g\) (grid/measurement height), \(u_*\), \(\theta_*\), prior flags (e.g., previous ON/OFF), time-of-day, season.

- **Anchored normalizations:**
  **Non-dimensional groups** (e.g., \(\zeta\)), **bounded scalings** (clip extreme gradients), **site-relative z-scaling** to mitigate \(\Delta z\) sensitivity.

- **Interaction terms:**
  **Shear–stratification ratio** \(\frac{S^2}{\Gamma}\), **bulk Richardson** proxies, **hysteresis features** (recent stability integral).

- **Feature hygiene:**
  **Collinearity pruning** (VIF checks), **monotonic sign enforcement** (e.g., ensure \(\Gamma>0\) for stable cases via absolute or max).

---

## Symbolic regression setup

- **Operator set (lean, physics-aware):**
  - **Add/multiply:** +, \(\times\)
  - **Divide with guard:** safe_div(x, y, eps)
  - **Unary transforms:** \(\log(1+x)\), \(\exp(x)\), \(|x|\)
  - **Bounded nonlinearity:** \(\tanh(x)\)
  - **Power (restricted):** \(x^\alpha\), \(\alpha \in [0.5, 2]\)

- **Complexity control:**
  **Description length penalty**, **max tree depth**, **max terms**, and **coefficient L1/L2** regularization.

- **Physics constraints:**
  **Monotonicity**: enforce/penalize desired signs w.r.t. \(\Gamma\) and \(S\); **output bounds**: hard clamp to \([0.15, 2.0]\).

- **Custom loss (multi-objective):**
  - **Regression term:** MAE on \(\mathrm{Ri_c^*}\).
  - **Event term:** 1 − F1 within \(\pm \Delta t\).
  - **Constraint penalties:** monotonicity violations, out-of-bounds before clamp, stiffness (large curvature) penalties.

- **Search regime:**
  **Two-stage**: coarse wide operator space → refined with pruned set; **Pareto tracking** of error vs complexity.

---

## Training protocol and selection

- **Cross-validation:**
  **By site/time blocks**; ensure regime diversity in train, hold-out in validation/test.

- **Curriculum:**
  **Start neutral/mild stability**, then include strong stable; **augment** with synthetic MOST to cover tails.

- **Pareto frontier:**
  **Select candidates** across complexity levels; prefer formulas with clear physical interpretation.

- **Checkpointing:**
  **Store expressions**, coefficients, metrics, and constraint diagnostics per candidate.

---

## Constraint enforcement and coefficient calibration

- **Clamping:**
  Apply \(\mathrm{Ri_c^*} \leftarrow \min(\max(\hat{f}(\mathbf{x}), 0.15), 2.0)\).

- **Monotonicity tests:**
  **Numerical partials**: finite differences on \(\Gamma\) and \(S\) grids; reject candidates violating required signs beyond tolerance.

- **Bounded regression:**
  **Refit coefficients** on held-out validation via constrained linear/nonlinear least squares within physical bounds.

- **Robustness to \(\Delta z\):**
  **Stress-test** with varying layer thickness; add penalty if sensitivity exceeds threshold.

---

## Validation and closed-loop tests

- **Statistical validation:**
  **MAE/RMSE** of \(\mathrm{Ri_c^*}\) and **F1/precision/recall** for event timing on held sites.

- **Closed-loop model test:**
  **Integrate** \(\mathrm{Ri_c^*}\) into the column scheme; check **stability**, **energy conservation**, **surface flux biases**, and **no crash** conditions across diurnal cycles and synoptic events.

- **Ablations:**
  **Feature drop** tests, **operator ablations**, and **baseline comparisons** against fixed \(\mathrm{Ri_c}\) and simple LUTs.

- **Sensitivity analyses:**
  **One-at-a-time** sweeps on \(\Gamma, S, \zeta, u_*\) to map response surfaces; confirm intended monotonicity and saturation.

---

## Model finalization and exports

- **Finalists:**
  **Pick simplest** expression meeting success criteria and safety checks; include **fallback** fixed \(\mathrm{Ri_c}\) path.

- **Runtime artifacts:**
  **Callable function** with guard rails, **LUT** for constrained regimes, and **ONNX** export if embedding in heterogeneous stacks.

- **Versioning:**
  **ModelVersion**, **DatasetVersion**, **Thresholds**, and **ConstraintConfig** embedded in metadata.

---

## Deployment recipe and logging

- **Inference speed:**
  **Benchmark** against budget; precompute reusable features (e.g., \(\zeta\)).

- **Integration points:**
  **Single entry function** `Ri_c_star(state)`; implement **safe_div**, **clamp**, **nan guards**, and **range checks**.

- **Telemetry:**
  **Structured logs** per call: inputs, \(\hat{f}\) pre-clamp, post-clamp, constraint flags, and runtime.

- **Fallback logic:**
  **Degrade gracefully** to fixed \(\mathrm{Ri_c}\) on invalid inputs, extreme regimes, or constraint violation detection.

---

## Reproducibility and artifacts

- **Environment:**
  **Containerize** (exact package versions), **seed RNGs**, and track **hardware**.

- **Run records:**
  **Hyperparameters**, operator sets, penalties, split indices, and data hashes.

- **Bundle:**
  **Formula text**, serialized function, **unit/integration tests**, **small reference dataset**, and **plots** (Pareto, response surfaces).

---

## Example snippets

#### Safe operators, clamp, and monotonicity checks

```python
import numpy as np

EPS = 1e-6
RI_MIN, RI_MAX = 0.15, 2.0

def safe_div(x, y, eps=EPS):
    return x / np.clip(y, -eps, eps)

def clamp(x, lo=RI_MIN, hi=RI_MAX):
    return np.minimum(np.maximum(x, lo), hi)

def ri_c_star(state, coeff):
    # Example symbolic form (candidate)
    Gamma = np.maximum(state["dtheta_dz"], 0.0)           # K/m
    S = np.maximum(state["dU_dz_abs"], 0.0)               # 1/s
    zeta = state["z_over_L"]                              # unitless
    dz = state["delta_z"]                                 # m

    term1 = coeff.a0 + coeff.a1 * np.log1p(Gamma) + coeff.a2 * np.log1p(S)
    term2 = coeff.a3 * np.tanh(coeff.a4 * zeta) + coeff.a5 * safe_div(Gamma, S + EPS)
    raw = term1 + term2 + coeff.a6 * np.tanh(coeff.a7 * dz / (1.0 + dz))
    return clamp(raw)

def monotonicity_violations(state_grid, coeff):
    # Finite-difference checks on Gamma and S
    vio = {"Gamma_pos": 0, "S_neg": 0}
    for Gamma in state_grid["Gamma"]:
        for S in state_grid["S"]:
            s = dict(dtheta_dz=Gamma, dU_dz_abs=S, z_over_L=state_grid["zeta"], delta_z=state_grid["dz"])
            base = ri_c_star(s, coeff)
            s["dtheta_dz"] = Gamma * 1.05
            dGamma = ri_c_star(s, coeff) - base
            s["dtheta_dz"] = Gamma
            s["dU_dz_abs"] = S * 1.05
            dS = ri_c_star(s, coeff) - base
            vio["Gamma_pos"] += int(dGamma < -1e-4)   # should not decrease with Gamma (stable)
            vio["S_neg"]    += int(dS     >  1e-4)   # should not increase with shear
    return vio
```

#### Custom multi-objective loss (pseudo-code)

```python
def loss(pred_Ric, true_Ric, pred_events, true_events, raw_unclamped, monotone_penalty_w):
    mae = np.mean(np.abs(pred_Ric - true_Ric))
    precision, recall = compute_precision_recall(pred_events, true_events)
    f1 = 2 * precision * recall / max(precision + recall, 1e-6)
    event_term = 1.0 - f1

    out_of_bounds = np.mean(np.maximum(0, RI_MIN - raw_unclamped) + np.maximum(0, raw_unclamped - RI_MAX))
    curvature = np.mean(np.abs(np.gradient(np.gradient(raw_unclamped))))
    monotone_pen = monotone_penalty_w * compute_monotone_penalty()

    return mae + event_term + 0.1 * out_of_bounds + 0.01 * curvature + monotone_pen
```

#### Deployment interface and logging

```python
class RiCModel:
    def __init__(self, coeff, meta):
        self.coeff = coeff
        self.meta = meta  # ModelVersion, DatasetVersion, Thresholds

    def infer(self, state):
        raw = self._raw(state)
        ric = clamp(raw)
        flags = self._constraints_check(state, raw, ric)
        self._log(state, raw, ric, flags)
        if flags["invalid_input"] or flags["hard_violation"]:
            return self.meta["fallback_Ri_c"]
        return ric

    # Implement _raw, _constraints_check, and _log with guards and telemetry
```

---

## Practical tips and pitfalls

- **Keep operators minimal:** Each extra operator bloats search space and risks overfitting; add only when payload is clear.
- **Test hysteresis explicitly:** Include prior ON/OFF flags to capture persistence but ensure bounded influence.
- **Bound gradients:** Penalize extreme curvature to avoid stiff coupling in the column solver.
- **Prefer saturation:** Use \(\tanh\) or \(\log(1+x)\) over unconstrained polynomials for stable behavior at tails.
- **Always include fallback:** A fixed \(\mathrm{Ri_c}\) path prevents degraded performance under unexpected regimes.

If you want, share your current feature list and event labeling logic, and I’ll tailor the operator set, loss weights, and monotonicity rules to your datasets and solver constraints.