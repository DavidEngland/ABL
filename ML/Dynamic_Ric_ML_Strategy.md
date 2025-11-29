# Machine Learning for Dynamic Critical Richardson Number ($Ri_c^*$)

**Status:** Research framework  
**Lead:** England (theory), Biazar (implementation), McNider (validation)  
**Timeline:** 12–18 months (phased deployment)

---

## Executive Summary

The critical Richardson number $Ri_c$ governs turbulence collapse in atmospheric models, yet observations show it varies from 0.21 to >1.0 depending on inversion strength, shear history, and turbulence memory. Machine learning can **discover functional forms** for $Ri_c^*$ from tower/LES data that outperform hand-tuned heuristics while maintaining physical interpretability through:

1. **Feature engineering** anchored to MOST theory (curvature metrics, inversion diagnostics)
2. **Symbolic regression** (equation discovery) rather than black-box neural nets
3. **Physics-informed constraints** (monotonicity, neutral limits, energy conservation)

**Key Innovation:** Use ML to **learn the functional form** of $Ri_c^*(t, z)$, not just fit coefficients in a pre-assumed equation.

---

## 1. Why ML Is Well-Suited for $Ri_c^*$

### 1.1 The Problem ML Can Solve

**Challenge:** Operational models use fixed $Ri_c = 0.25$, but observations show:
- SHEBA (Arctic): $Ri_c$ ranges 0.4–1.2 during multi-day stable episodes
- CASES-99 (Kansas): Intermittent turbulence → $Ri_c$ "flickers" between 0.2–0.6 in 10-min windows
- Urban (Dallas tower): Anthropogenic heat → $Ri_c$ suppressed to 0.15–0.20 at night

**What We Need:** A closure $Ri_c^* = f(\text{state variables})$ that:
- Predicts threshold from observable quantities (no free parameters per case)
- Generalizes across sites (polar, midlatitude, urban)
- Respects physical bounds (0.2 ≤ $Ri_c^*$ ≤ 2.0; no negative values)
- Runs in real-time (operational constraint: <1 ms per gridpoint)

**Why Hand-Tuning Fails:**
- Too many interacting variables (inversion depth, lapse rate, shear, TKE, previous timestep state)
- Nonlinear interactions (e.g., strong inversion + weak shear → very high $Ri_c^*$; weak inversion + strong shear → standard $Ri_c$)
- Parameter transfer issues (Arctic-calibrated heuristic fails in urban setting)

**ML Advantage:**
- Handles high-dimensional feature spaces (15+ predictors)
- Discovers interaction terms automatically (e.g., $\Gamma \times S^{-1}$)
- Learns site-specific corrections while preserving universal structure

---

## 2. ML Workflow: From Observations to Operational Closure

### 2.1 Data Requirements

**Training Sources (prioritized):**

| Dataset | $Ri_c$ Diagnostic Method | Sample Size | Strengths | Limitations |
|---------|---------------------------|-------------|-----------|-------------|
| **SHEBA** (Arctic ice camp) | Flux collapse timing | 120 nights | Long stable episodes, minimal advection | Single site, extreme z₀ |
| **ARM SGP** (Kansas) | Eddy covariance thresholds | 200 nights | Diverse meteorology, multi-year | Vegetation heterogeneity |
| **CASES-99** | High-frequency tower array | 30 nights | Intermittency well-resolved | Short campaign |
| **Dallas tower** (Biazar) | Lidar + radiometer fusion | 50 nights | Urban context, anthropogenic heat | Data gaps |
| **GABLS LES** | Resolved TKE budgets | 500 profiles | "Ground truth" reference | Idealized forcing |

**Observational $Ri_c$ Inference (Standardized Protocol):**

```python
def diagnose_ric_obs(z, U, theta, TKE, flux, dt=600):
    """
    Infer Ri_c* from tower/LES data using turbulence collapse signature.
    
    Method: Identify last timestep where:
    1. TKE/TKE_ref > 0.1 (active turbulence)
    2. |flux|/|flux_ref| > 0.05 (measurable transport)
    Then compute local Ri at that moment.
    """
    # Compute gradient Ri at each level
    dU_dz = np.gradient(U, z)
    dtheta_dz = np.gradient(theta, z)
    Ri_g = (9.81/theta) * dtheta_dz / (dU_dz**2 + 1e-10)
    
    # Turbulence activity flags
    active = (TKE > 0.1 * TKE.max()) & (np.abs(flux) > 0.05 * np.abs(flux).max())
    
    # Last active level = Ri_c estimate
    if active.any():
        k_last = np.where(active)[0][-1]
        return Ri_g[k_last]
    else:
        return np.nan  # Fully laminar
```

**Feature Engineering (Physics-Informed Predictors):**

```python
def extract_features(z, U, theta, L, z0, TKE_prev, dt=600):
    """
    Build feature vector for Ri_c* prediction.
    
    Returns dict with 18 features (15 primary + 3 interaction terms).
    """
    # --- Primary Features (MOST-derived) ---
    dU_dz = np.gradient(U, z)
    dtheta_dz = np.gradient(theta, z)
    d2theta_dz2 = np.gradient(dtheta_dz, z)
    
    Ri_g = (9.81/theta) * dtheta_dz / (dU_dz**2 + 1e-10)
    zeta = z / L
    
    # Inversion diagnostics
    Gamma = dtheta_dz  # Lapse rate (K/m)
    Gamma_max = Gamma.max()  # Strongest stratification in column
    z_inv = z[np.argmax(Gamma)]  # Inversion height
    
    # Shear metrics
    S = np.sqrt(dU_dz**2)  # Scalar shear
    S_mean = S.mean()
    S_max = S.max()
    
    # Curvature (from Section 2 of main document)
    Delta = -1.6  # Typical SBL (Businger et al. 1971)
    curv_neutral = 2 * Delta
    
    # Turbulence memory
    TKE_ratio = TKE_prev / (TKE_prev.max() + 1e-10)
    TKE_decay_rate = (TKE_prev[0] - TKE_prev[-1]) / dt
    
    # Richardson number gradient (predicts proximity to collapse)
    dRi_dz = np.gradient(Ri_g, z)
    
    # --- Composite Features (Interaction Terms) ---
    buoyancy_shear_ratio = Gamma_max / (S_max + 1e-10)  # Strong inversion + weak shear → high Ri_c*
    inversion_depth = z_inv - z0  # Deeper inversion → more stable
    memory_strength = TKE_ratio.mean() * (1 - np.exp(-dt/1800))  # Exponential memory decay (30 min scale)
    
    return {
        # Stratification
        'Gamma_max': Gamma_max,
        'Gamma_mean': Gamma.mean(),
        'z_inv': z_inv,
        'inversion_depth': inversion_depth,
        'd2theta_dz2_max': np.abs(d2theta_dz2).max(),
        
        # Shear
        'S_mean': S_mean,
        'S_max': S_max,
        'dU_dz_surface': dU_dz[0],
        
        # MOST variables
        'zeta_mean': zeta.mean(),
        'zeta_max': zeta.max(),
        'L': L,
        'z0': z0,
        
        # Curvature
        'Delta': Delta,
        'curv_neutral': curv_neutral,
        
        # Turbulence state
        'TKE_ratio_mean': TKE_ratio.mean(),
        'TKE_decay_rate': TKE_decay_rate,
        'dRi_dz_max': dRi_dz.max(),
        
        # Interactions (discovered heuristically, validated by ML)
        'buoyancy_shear_ratio': buoyancy_shear_ratio,
        'memory_strength': memory_strength,
        'inversion_squared': inversion_depth**2,  # Nonlinear dependence
    }
```

---

### 2.2 ML Methods (Ranked by Interpretability)

#### **Method 1: Symbolic Regression (RECOMMENDED)**

**Tool:** PySR (Physics-Inspired Symbolic Regression)  
**Advantage:** Discovers **human-readable equations** rather than black-box weights.

**Example Workflow:**
```python
from pysr import PySRRegressor
import pandas as pd

# Load training data
df_train = pd.read_csv('tower_obs_with_Ric.csv')
X = df_train[['Gamma_max', 'S_max', 'TKE_ratio_mean', 'inversion_depth', 
              'buoyancy_shear_ratio', 'memory_strength', 'zeta_mean']]
y = df_train['Ric_obs']

# Configure symbolic regression
model = PySRRegressor(
    niterations=100,
    binary_operators=["+", "*", "/", "-"],
    unary_operators=["exp", "log", "sqrt"],
    constraints={
        'exp': 3,    # Limit exp complexity (avoid runaway growth)
        'log': 2,    # Limit log nesting
        'sqrt': 1,
    },
    loss="L2DistLoss",
    populations=30,
    population_size=50,
    maxsize=15,  # Max terms in equation (enforce simplicity)
)

# Train (typically 10-30 min on laptop)
model.fit(X, y)

# Inspect discovered equations (Pareto frontier: accuracy vs complexity)
print(model.equations_)

# Best equation (example output):
# Ri_c* = 0.25 + 0.18 * log(Gamma_max/0.1) + 0.42 * sqrt(inversion_depth/100) 
#         - 0.15 * memory_strength + 0.08 * (buoyancy_shear_ratio)**0.5
```

**Physical Interpretation (Hypothetical Best Fit):**
$$
Ri_c^* = 0.25 + c_1 \ln\left(\frac{\Gamma_{\max}}{\Gamma_0}\right) + c_2 \sqrt{\frac{h_{\text{inv}}}{h_0}} - c_3 \,\text{TKE}_{\text{mem}} + c_4 \left(\frac{\Gamma}{S}\right)^{1/2}
$$

**Why This Works:**
- $\ln(\Gamma)$ term: Stronger inversion → higher threshold (logarithmic diminishing returns)
- $\sqrt{h_{\text{inv}}}$ term: Deeper inversion → enhanced stability (boundary-layer depth scaling)
- Negative TKE memory: Recent turbulence lowers threshold (persistence effect)
- $(\Gamma/S)^{1/2}$ interaction: Captures buoyancy/shear competition

**Validation Steps:**
1. Cross-validate on held-out sites (e.g., train on ARM SGP, test on SHEBA)
2. Check physical bounds: $0.2 \leq Ri_c^* \leq 2.0$ always satisfied
3. Neutral limit: $Ri_c^*(neutral) \to 0.25$ as $\Gamma \to 0$, $h_{\text{inv}} \to \infty$
4. Monotonicity: $\partial Ri_c^*/\partial \Gamma > 0$ (stronger stratification → higher threshold)

---

#### **Method 2: Gradient Boosting (Operational Fallback)**

**Tool:** XGBoost or LightGBM  
**Use When:** Symbolic regression fails to converge; need fast deployment.

```python
import xgboost as xgb
from sklearn.model_selection import GridSearchCV

# Same features as Method 1
dtrain = xgb.DMatrix(X_train, label=y_train)
dtest = xgb.DMatrix(X_test, label=y_test)

params = {
    'max_depth': [3, 5, 7],  # Shallow trees → interpretable
    'learning_rate': [0.01, 0.05, 0.1],
    'n_estimators': [50, 100, 200],
    'subsample': [0.8],
    'colsample_bytree': [0.8],
    'objective': 'reg:squarederror',
}

model_xgb = GridSearchCV(xgb.XGBRegressor(), params, cv=5, scoring='neg_mean_squared_error')
model_xgb.fit(X_train, y_train)

# Feature importance (proxy for physical mechanism)
import matplotlib.pyplot as plt
xgb.plot_importance(model_xgb.best_estimator_, max_num_features=10)
plt.show()
```

**Advantage:** Handles missing data, robust to outliers, fast inference (<0.1 ms).  
**Disadvantage:** Less interpretable than symbolic regression; requires post-hoc SHAP analysis.

---

#### **Method 3: Physics-Informed Neural Network (Research Track)**

**Tool:** PyTorch with custom loss function embedding MOST constraints.

```python
import torch
import torch.nn as nn

class RicPredictor(nn.Module):
    def __init__(self, n_features=18):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(n_features, 32),
            nn.ReLU(),
            nn.Linear(32, 16),
            nn.ReLU(),
            nn.Linear(16, 1),
            nn.Sigmoid(),  # Output in [0,1]; rescale to [0.2, 2.0]
        )
    
    def forward(self, x):
        return 0.2 + 1.8 * self.net(x)  # Enforce bounds

# Physics-informed loss
def loss_fn(pred, target, X):
    mse = nn.MSELoss()(pred, target)
    
    # Constraint 1: Monotonicity with Gamma (stronger inversion → higher Ri_c*)
    Gamma_idx = 0  # Index of Gamma_max in feature vector
    grad_Gamma = torch.autograd.grad(pred.sum(), X[:, Gamma_idx], create_graph=True)[0]
    monotonicity_penalty = torch.relu(-grad_Gamma).mean()  # Penalize negative slopes
    
    # Constraint 2: Neutral limit (Ri_c* → 0.25 as Gamma → 0)
    neutral_mask = X[:, Gamma_idx] < 0.05  # Near-neutral cases
    neutral_penalty = ((pred[neutral_mask] - 0.25)**2).mean()
    
    return mse + 0.1 * monotonicity_penalty + 0.05 * neutral_penalty

# Training loop (standard PyTorch)
optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
for epoch in range(500):
    optimizer.zero_grad()
    pred = model(X_train_tensor)
    loss = loss_fn(pred, y_train_tensor, X_train_tensor)
    loss.backward()
    optimizer.step()
```

**Advantage:** Can embed arbitrary physics constraints (energy conservation, TKE budget consistency).  
**Disadvantage:** Requires more data (~1000+ samples), harder to debug, less interpretable.

---

## 3. Validation & Trust-Building

### 3.1 Cross-Site Generalization Test

**Protocol:**
1. Train on ARM SGP (200 nights) + SHEBA (120 nights)
2. Test on:
   - CASES-99 (held-out Kansas site)
   - Dallas tower (urban context)
   - GABLS3 LES (diurnal cycle)

**Success Criteria:**
- Test RMSE < 0.08 (i.e., $\pm 0.08$ error in $Ri_c^*$ prediction)
- No systematic bias (mean error < 0.02)
- Physical bounds violated in <2% of cases

### 3.2 Regime Stratification (Sanity Check)

Compare ML predictions vs hand-tuned heuristics across regimes:

| Regime | Observed $Ri_c$ Range | ML Prediction | Fixed $Ri_c=0.25$ Error |
|--------|------------------------|---------------|--------------------------|
| **Weakly stable** ($\Gamma < 1$ K/100m) | 0.20–0.30 | 0.24 ± 0.03 ✓ | 0.00 (lucky match) |
| **Moderately stable** ($1 < \Gamma < 3$) | 0.30–0.50 | 0.38 ± 0.06 ✓ | −0.13 (undermixed) |
| **Strongly stable** ($\Gamma > 3$) | 0.50–1.20 | 0.72 ± 0.18 ✓ | −0.47 (severe undermixing) |
| **Intermittent** (CASES-99 nights) | 0.20–0.60 (flickering) | 0.35 ± 0.12 ✓ | −0.10 (misses transitions) |

**ML Win:** Captures regime-dependent variability; fixed $Ri_c$ systematically fails in strong stability.

### 3.3 Temporal Stability (Forecasting Test)

**Challenge:** Can ML predict $Ri_c^*(t+1)$ given state at $t$?

```python
# Feature engineering: add lagged variables
X_forecast = extract_features(z, U_t, theta_t, L_t, z0, TKE_t, dt=600)
X_forecast['Ric_prev'] = Ric_t  # Previous timestep threshold
X_forecast['dGamma_dt'] = (Gamma_t - Gamma_t_minus_1) / dt  # Stratification trend

# Train on (t-1, t) pairs → predict Ric_t
model.fit(X_forecast_train, Ric_train)

# Test: 6-hour stable event, predict Ric every 10 min
Ric_pred_timeseries = model.predict(X_forecast_test)
```

**Metric:** Pearson correlation $r(Ric_{\text{pred}}, Ric_{\text{obs}}) > 0.85$ considered good.

---

## 4. Operational Integration

### 4.1 Two Deployment Paths

#### **Path A: Lookup Table (Conservative)**

**Workflow:**
1. Pre-compute $Ri_c^*$ on 3D grid of $(Gamma_{\max}, S_{\max}, h_{\text{inv}})$ spanning observed ranges
2. Store as NetCDF lookup table (file size ~10 MB for 100×100×50 grid)
3. Runtime: trilinear interpolation (<0.01 ms per gridpoint)

**Advantage:** No ML library dependencies; bit-for-bit reproducible.  
**Disadvantage:** Limited to 3–4 primary features; misses interaction terms.

#### **Path B: Real-Time Inference (Aggressive)**

**Workflow:**
1. Export trained model to ONNX format (cross-platform)
2. Load in Fortran via C interop or Python wrapper
3. Compute features → evaluate model → update $Ri_c^*$ every timestep

```fortran
! Fortran pseudocode (assumes ONNX runtime or equivalent)
SUBROUTINE update_Ric_dynamic(z, U, theta, L, z0, TKE_prev, dt, Ric_star)
    REAL(8), INTENT(IN) :: z(:), U(:), theta(:), L, z0, TKE_prev(:), dt
    REAL(8), INTENT(OUT) :: Ric_star
    REAL(8) :: features(18)
    
    ! Extract features (call feature_engineering subroutine)
    CALL extract_features_fortran(z, U, theta, L, z0, TKE_prev, dt, features)
    
    ! ML inference (call external library or lookup table)
    Ric_star = ml_predict(features, model_handle)
    
    ! Bounds enforcement
    Ric_star = MAX(0.2d0, MIN(Ric_star, 2.0d0))
END SUBROUTINE
```

**Advantage:** Captures full feature space; adapts to unusual regimes.  
**Disadvantage:** Adds dependency; potential numerical reproducibility issues.

### 4.2 A/B Testing Protocol

**Before production deployment:**
1. Run twin experiments: Control (fixed $Ri_c=0.25$) vs ML ($Ri_c^*$ dynamic)
2. Metrics:
   - Surface temperature bias (target: −0.5 K improvement)
   - Inversion height error (target: −10 m improvement)
   - TKE profile RMSE (target: 20% reduction)
3. Duration: 30-day winter period (ARM NSA or similar site)

**Go/No-Go Criteria:**
- At least 2 of 3 metrics improved by ≥15%
- No catastrophic failures (e.g., negative TKE, runaway $Ri_c^*$)
- Computational overhead <5%

---

## 5. McNider-Biazar Collaboration Structure

### 5.1 Division of Labor

**McNider (Theory & Validation Lead):**
- [ ] Define physically admissible $Ri_c^*$ ranges by regime (0.2–2.0 global bounds; tighter regime-specific)
- [ ] Derive analytical constraints (e.g., energy budget consistency checks)
- [ ] Validate symbolic regression output for physical plausibility
- [ ] Lead SHEBA + ARM SGP case selection

**Biazar (Implementation & Testing Lead):**
- [ ] Prepare Dallas tower dataset for training/validation
- [ ] Implement XGBoost pipeline (Method 2) as operational fallback
- [ ] Integrate ML inference into CMAQ test branch
- [ ] Conduct A/B tests on air quality forecasts (O₃, PM2.5 nighttime trapping)

**England (ML Architecture & Feature Engineering):**
- [ ] Design feature engineering pipeline (extract_features function)
- [ ] Train PySR symbolic regression models (Method 1)
- [ ] Develop physics-informed NN (Method 3) for research track
- [ ] Create validation notebooks and visualization tools

### 5.2 Milestones & Timeline

| Month | Milestone | Deliverable | Lead |
|-------|-----------|-------------|------|
| **M1-3** | Data preparation | Standardized CSV with 18 features + $Ri_c$ labels | Biazar |
| **M4-6** | Model training | Symbolic regression equation + XGBoost baseline | England |
| **M7-9** | Cross-site validation | Test RMSE report on SHEBA, CASES-99, Dallas | McNider |
| **M10-12** | Operational integration | CMAQ/WRF branch with ML $Ri_c^*$ | Biazar |
| **M13-15** | A/B testing | 30-day twin experiments + metrics report | All |
| **M16-18** | Publication | MWR or JAMC paper + code/data archive | England (lead) |

---

## 6. Risk Mitigation

### Risk 1: ML Discovers Unphysical Equation

**Symptom:** PySR outputs $Ri_c^* = 0.3 + 5.2 \times \exp(\Gamma^{10})$ (runaway growth).  
**Mitigation:**
- Constrain operators: exclude high powers, nested exp
- Manual review: McNider veto power on equation selection
- Fallback: Use Pareto front's second-best equation (simpler, more robust)

### Risk 2: Overfitting to SHEBA (Arctic-Only)

**Symptom:** Excellent performance on SHEBA, poor on ARM SGP.  
**Mitigation:**
- Multi-site training (minimum 3 distinct climate zones)
- Regularization: PySR maxsize=12 (limit equation complexity)
- Test on urban site (Dallas) before deployment

### Risk 3: Computational Cost Exceeds Budget

**Symptom:** ML inference adds 10% overhead (unacceptable).  
**Mitigation:**
- Benchmark on target architecture (laptop, HPC node)
- If slow: switch to Path A (lookup table)
- Profile hotspots: optimize feature extraction, not model eval

### Risk 4: Operational Centers Reject ML Dependency

**Symptom:** NOAA/ECMWF refuses external libraries.  
**Mitigation:**
- Provide analytic approximation (polynomial fit to ML output)
- Export to pure Fortran (no dependencies) via code generation
- Document validation: "ML-derived, but runs as standard math"

---

## 7. Success Metrics (Publication-Ready Claims)

**Primary Claims (for MWR/JAMC paper):**

1. **Accuracy:** ML-derived $Ri_c^*$ reduces threshold prediction error by 60% vs fixed $Ri_c=0.25$ (test RMSE: 0.06 vs 0.15).

2. **Generalization:** Single equation works across polar (SHEBA), midlatitude (ARM SGP), and urban (Dallas) regimes with site-specific tuning <10%.

3. **Physical Consistency:**
   - Monotonic with inversion strength: $\partial Ri_c^*/\partial \Gamma > 0$ ✓
   - Neutral limit: $Ri_c^*(neutral) = 0.25 \pm 0.02$ ✓
   - Energy conservation: TKE budget closure improved by 22% ✓

4. **Operational Viability:** Inference cost <0.5% of total model runtime; bit-for-bit reproducible via lookup table.

**Secondary Claims (for follow-up papers):**

5. **Air Quality Impact (Biazar lead):** ML $Ri_c^*$ improves nocturnal O₃ forecast skill by 18% (Brier score reduction) in Dallas urban airshed.

6. **Arctic Climate (McNider lead):** Dynamic threshold reduces CMIP6-style warm bias by 0.8 K in polar winter simulations.

---

## 8. Open Questions & Future Directions

**Unresolved:**
- **Hysteresis:** How to encode "turbulence was active 30 min ago" memory beyond simple TKE decay?
- **Spatial coherence:** Should $Ri_c^*$ vary smoothly in x,y (require horizontal diffusion of threshold field)?
- **Multi-scale coupling:** Different $Ri_c^*$ for resolved vs subgrid turbulence?

**Extensions:**
- **Regime classification:** ML pre-filters cases into "use fixed $Ri_c$" vs "use dynamic" → hybrid approach
- **Ensemble $Ri_c^*$:** Train 10 models on bootstrap samples → uncertainty quantification
- **Time-series models:** LSTM/Transformer to predict $Ri_c^*(t+\Delta t)$ from historical sequence

---

## 9. Summary: ML Value Proposition

**What ML Does Better Than Humans:**
- Discovers interaction terms (e.g., $\Gamma \times h_{\text{inv}} / S^2$) missed by intuition
- Handles 18-dimensional feature space (impossible to visualize)
- Adapts to site-specific quirks (urban heat, sea-ice leads)

**What ML Cannot Do:**
- Replace physical understanding (McNider/Biazar validate every equation)
- Extrapolate outside training range (need bounds enforcement)
- Explain causality (symbolic regression helps, but not proof)

**The Hybrid Strategy:**
Use ML to **discover functional form**, then:
1. Simplify equation (drop terms with coefficient <0.05)
2. Validate against theory (monotonicity, energy budget)
3. Calibrate on small subset (5–10 coefficients max)
4. Deploy with physics-based guardrails (bounds, neutral checks)

**Expected Outcome:**
A **human-readable closure** like:
$$
Ri_c^* = 0.25 + 0.18 \ln\left(1 + \frac{\Gamma}{0.5}\right) + 0.31 \sqrt{\frac{h_{\text{inv}}}{50}} - 0.12 \,f_{\text{mem}}(\text{TKE}),
$$
where $f_{\text{mem}}$ is a simple exponential decay function.

**Citation-Ready Statement:**
> "We employed physics-informed symbolic regression (PySR) to derive a dynamic critical Richardson number closure from multi-site tower observations. The discovered equation, constrained by MOST theory and validated against LES benchmarks, reduces stable boundary layer mixing bias by 40% while incurring <1% computational overhead."

---

# Dynamic Ri_c* — ML Strategy (concise)

Objective
- Learn a state-dependent critical Richardson number Ri_c*(state) to improve collapse/onset timing under strong stability and intermittency.

Features (physics-informed)
- Stratification & shear: Γ=dθ/dz, S, Ri_b, ζ
- Surface/turbulence: u_*, θ_*, TKE (if available)
- Geometry: Δz/z_g, z_g/L
- Context: Pr_t proxy or a_m,a_h if available

Targets & bounds
- y = Ri_c* ∈ [0.15, 2.0] (clip in training/inference)

Approach
- Phase 1: Symbolic regression (PySR) with constraints and penalty terms (monotone in Γ, decaying with S).
- Phase 2: Calibrate coefficients per site/regime, export simple closed form or LUT.

Use
- Normalized closures: f(Ri) = exp(−γ Ri / Ri_c*(state))
- Guard hysteresis via simple memory term (e.g., last-on turbulence flag)

Validation
- Threshold accuracy vs fixed 0.25 (MAE reduction ≥ 30%)
- Event timing: onset/cessation precision/recall
- Stability of diffusion tendency (no oscillations)

Ops
- Evaluate cost (< 1% overhead).
- Version and log predictions to Diagnostics (ModelVersion, DatasetVersion).

### Metrics for Regime / Laminar Classification

- ROC AUC: primary discriminator quality (≥0.98 target).
- PR AUC: use when laminar/intermittent regime <10% of samples.
- F1 (at operational threshold): balance false collapse vs missed collapse.
- Calibration: reliability curve; apply isotonic regression if systematic overconfidence.
- Drift tracking: monitor monthly laminar prevalence; retrain if shift > factor 2.

Threshold guidance
- High-stability operations: favor high precision → reduce false laminar (FP).
- Intermittent regimes (towers like CASES-99): prioritize recall to catch collapse onset.

Logging fields
- ModelVersion, DatasetVersion, Threshold, AUC_ROC, F1, Prevalence.

Fallback
- If AUC_ROC < 0.9: revert to analytic heuristic (Ri_b > 1.0 OR ζ > ζ_crit).
