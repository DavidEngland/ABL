# Machine Learning for Richardson Number Curvature Corrections

**Purpose:** Augment physics-based MOST corrections with lightweight ML to handle:
- Fast bias estimation on coarse grids
- Adaptive correction parameter tuning
- Automatic detection of problematic stability regimes
- Real-time quality control for tower/model data

**Philosophy:** ML as a **diagnostic assistant**, not a physics replacement. All training uses MOST-consistent synthetic data or validated observations.

---

## 1. Core ML Applications (Ordered by Complexity)

### 1.1 ML Bias Predictor (Simplest — Start Here)

**Problem:** Computing the exact bias ratio $B = Ri_g(z_g) / Ri_b$ requires analytic integration of $\phi$ functions, which may be expensive or unavailable at runtime.

**ML Solution:** Train a regression model to predict $B$ directly from easily-available layer diagnostics.

**Inputs (Features):**
```python
features = [
    'delta_z',        # Layer thickness (m)
    'z_bottom',       # Bottom height (m)
    'z0',             # Surface roughness (m)
    'Ri_b_raw',       # Bulk Richardson (raw, from model)
    'zeta_geom',      # ζ at geometric mean height
    'dtheta_dz',      # Temperature gradient (K/m)
    'du_dz'           # Wind shear (1/s)
]
```

**Target (Label):**
```python
target = 'B_true'  # Exact bias ratio from analytic MOST
```

**Training Data Generation (Synthetic):**
```python
import numpy as np
from scipy.integrate import quad

def generate_training_data(n_samples=10000):
    """Generate synthetic profiles with known B."""
    data = []
    
    # Parameter ranges (stable conditions only)
    delta_z_range = np.logspace(0.5, 2.5, 50)  # 3-300 m
    L_range = np.logspace(0.5, 2.5, 50)        # 3-300 m
    z0_range = np.logspace(-3, 0, 20)          # 0.001-1 m
    
    # MOST parameters (linear stable)
    a_m, a_h = 4.7, 7.8  # Businger-Dyer SBL
    
    for _ in range(n_samples):
        z0 = np.random.choice(z0_range)
        z_bottom = np.random.uniform(z0 + 1, 50)
        delta_z = np.random.choice(delta_z_range)
        z_top = z_bottom + delta_z
        z_geom = np.sqrt(z_bottom * z_top)
        
        L = np.random.choice(L_range)
        zeta_geom = z_geom / L
        
        # Skip if outside valid MOST range
        if zeta_geom > 2.0 or zeta_geom < 0.01:
            continue
        
        # Analytic Ri_g at geometric mean
        phi_m = 1 + a_m * zeta_geom
        phi_h = 1 + a_h * zeta_geom
        Ri_g_geom = zeta_geom * phi_h / phi_m**2
        
        # Analytic Ri_b (integrate Ri_g over layer)
        def integrand(z):
            zeta = z / L
            pm = 1 + a_m * zeta
            ph = 1 + a_h * zeta
            return zeta * ph / pm**2
        
        Ri_b_true, _ = quad(integrand, z_bottom, z_top)
        Ri_b_true /= delta_z
        
        # Bias ratio (target)
        B_true = Ri_g_geom / Ri_b_true if Ri_b_true > 1e-6 else 1.0
        
        # Synthetic bulk Ri (what model would compute)
        theta_bottom = 280.0
        theta_top = theta_bottom + 5.0 * (zeta_geom**0.5)  # Simplified
        dtheta = theta_top - theta_bottom
        
        U_bottom = 5.0
        u_star = 0.3
        kappa = 0.4
        U_top = U_bottom + (u_star / kappa) * np.log(z_top / z_bottom)
        dU = U_top - U_bottom
        
        g = 9.81
        theta_ref = 0.5 * (theta_bottom + theta_top)
        Ri_b_raw = (g / theta_ref) * (dtheta * delta_z) / (dU**2)
        
        dtheta_dz = dtheta / delta_z
        du_dz = dU / delta_z
        
        data.append({
            'delta_z': delta_z,
            'z_bottom': z_bottom,
            'z0': z0,
            'Ri_b_raw': Ri_b_raw,
            'zeta_geom': zeta_geom,
            'dtheta_dz': dtheta_dz,
            'du_dz': du_dz,
            'B_true': B_true
        })
    
    return pd.DataFrame(data)
```

**Model Training (Random Forest — Robust & Interpretable):**
```python
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_absolute_error, r2_score
import pandas as pd

# Generate data
df = generate_training_data(n_samples=10000)

# Features and target
X = df[['delta_z', 'z_bottom', 'z0', 'Ri_b_raw', 'zeta_geom', 
        'dtheta_dz', 'du_dz']]
y = df['B_true']

# Split
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=42
)

# Train
model_B = RandomForestRegressor(
    n_estimators=100,
    max_depth=10,
    min_samples_leaf=5,
    random_state=42
)
model_B.fit(X_train, y_train)

# Validate
y_pred = model_B.predict(X_test)
print(f"MAE: {mean_absolute_error(y_test, y_pred):.4f}")
print(f"R²: {r2_score(y_test, y_pred):.4f}")

# Feature importance
importances = pd.DataFrame({
    'feature': X.columns,
    'importance': model_B.feature_importances_
}).sort_values('importance', ascending=False)
print(importances)
```

**Deployment:**
```python
def predict_bias_ratio(model, delta_z, z_bottom, z0, Ri_b_raw, 
                       zeta_geom, dtheta_dz, du_dz):
    """Fast ML-based bias prediction."""
    features = np.array([[delta_z, z_bottom, z0, Ri_b_raw, 
                          zeta_geom, dtheta_dz, du_dz]])
    B_ml = model.predict(features)[0]
    return np.clip(B_ml, 1.0, 3.0)  # Physical bounds
```

---

### 1.2 Automated Correction Parameter Tuning (α, γ)

**Problem:** The McNider correction factor $f_c$ depends on tuning parameters $\alpha$ and $\gamma$ that should vary with stability regime and model resolution.

**ML Solution:** Learn optimal $(\alpha, \gamma)$ from validation against LES or tower observations.

**Training Setup:**
- **Inputs:** Grid properties, stability state, inversion strength
- **Targets:** $\alpha_{\text{opt}}$, $\gamma_{\text{opt}}$ that minimize flux RMSE
- **Training data:** Paired coarse-grid runs vs fine-grid/LES truth

**Objective Function (Offline Optimization):**
```python
def optimize_alpha_gamma(coarse_model_output, fine_truth, 
                         stability_regime):
    """Find (α, γ) minimizing surface flux error."""
    from scipy.optimize import minimize
    
    def loss(params):
        alpha, gamma = params
        # Apply correction with these parameters
        K_corrected = apply_mcnider_correction(
            coarse_model_output['K'], 
            coarse_model_output['zeta'],
            coarse_model_output['delta_z'],
            alpha=alpha, gamma=gamma
        )
        # Compute flux and compare
        flux_corrected = compute_surface_flux(K_corrected)
        flux_true = fine_truth['surface_flux']
        return np.mean((flux_corrected - flux_true)**2)
    
    result = minimize(loss, x0=[1.0, 1.5], 
                     bounds=[(0.5, 2.5), (0.5, 2.5)])
    return result.x
```

**ML Meta-Model (Learn the Optimizer):**
```python
# Features: regime characteristics
features_regime = [
    'mean_zeta',           # Average stability
    'max_zeta',            # Peak stability
    'inversion_strength',  # Δθ across inversion
    'delta_z_mean',        # Mean layer thickness
    'z0'                   # Surface roughness
]

# Targets: optimal parameters
targets = ['alpha_opt', 'gamma_opt']

# Train on many synthetic runs
model_params = RandomForestRegressor(n_estimators=50)
model_params.fit(X_regime_train, y_params_train)

# Deploy: instant parameter prediction
alpha_ml, gamma_ml = model_params.predict([[mean_zeta, max_zeta, 
                                            inv_str, dz_mean, z0]])[0]
```

---

### 1.3 Problematic Layer Detection (Classification)

**Problem:** Not all layers need correction — applying $f_c$ everywhere wastes computation and can introduce noise in neutral/unstable regions.

**ML Solution:** Binary classifier to flag layers where curvature bias is significant.

**Training Labels:**
```python
def label_layer(B, zeta, Ri_b):
    """Rule-based labeling for training."""
    if B > 1.15 and zeta > 0.1 and Ri_b > 0.05:
        return 1  # Apply correction
    else:
        return 0  # Skip correction
```

**Model (Logistic Regression for Interpretability):**
```python
from sklearn.linear_model import LogisticRegression

# Features
X_class = df[['B_true', 'zeta_geom', 'Ri_b_raw', 'delta_z']]
y_class = df.apply(lambda row: label_layer(row['B_true'], 
                                           row['zeta_geom'], 
                                           row['Ri_b_raw']), axis=1)

clf = LogisticRegression(penalty='l2', C=1.0)
clf.fit(X_class, y_class)

# Interpretation: coefficients show which features matter most
print("Coefficients:", dict(zip(X_class.columns, clf.coef_[0])))
```

**Runtime Use:**
```python
def should_apply_correction(clf, B_pred, zeta, Ri_b, delta_z):
    """Fast decision: correct this layer?"""
    prob = clf.predict_proba([[B_pred, zeta, Ri_b, delta_z]])[0, 1]
    return prob > 0.7  # High confidence threshold
```

---

### 1.4 Quality Control for Tower Data

**Problem:** Real observations have gaps, outliers, and instrument errors that corrupt $Ri_g$ calculations.

**ML Solution:** Anomaly detection + imputation for $\partial\theta/\partial z$ and $\partial U/\partial z$.

**Approach A: Isolation Forest (Outlier Detection):**
```python
from sklearn.ensemble import IsolationForest

# Train on "clean" synthetic profiles
iso_forest = IsolationForest(contamination=0.05, random_state=42)
iso_forest.fit(X_clean)  # Features: dtheta_dz, du_dz, Ri_g, etc.

# Flag outliers in tower data
outliers = iso_forest.predict(X_tower) == -1
print(f"Flagged {outliers.sum()} of {len(X_tower)} observations")
```

**Approach B: Gradient Imputation (Gaussian Process Regression):**
```python
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel

# Impute missing dtheta/dz using nearby levels and time context
kernel = ConstantKernel(1.0) * RBF([10.0, 1.0])  # Spatial + temporal
gp = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=5)

# Train on valid observations
X_valid = tower_data[~tower_data['dtheta_dz'].isna()][['z', 'time_hour']]
y_valid = tower_data.loc[~tower_data['dtheta_dz'].isna(), 'dtheta_dz']
gp.fit(X_valid, y_valid)

# Predict missing values
X_missing = tower_data[tower_data['dtheta_dz'].isna()][['z', 'time_hour']]
y_imputed = gp.predict(X_missing)
```

---

## 2. Advanced Applications

### 2.1 Emulating Fine-Grid Behavior

**Goal:** Given a coarse-grid profile, predict what the fine-grid $Ri_g(z)$ would have looked like.

**Architecture:** Neural network (simple feedforward, 2-3 hidden layers).

**Training:**
```python
import torch
import torch.nn as nn

class RiProfileEmulator(nn.Module):
    def __init__(self, n_features=7, n_layers=32):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(n_features, n_layers),
            nn.ReLU(),
            nn.Linear(n_layers, n_layers),
            nn.ReLU(),
            nn.Linear(n_layers, n_layers // 2),
            nn.ReLU(),
            nn.Linear(n_layers // 2, n_layers)  # Output: fine-grid Ri values
        )
    
    def forward(self, x):
        return self.net(x)

# Training loop (simplified)
model = RiProfileEmulator()
optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
criterion = nn.MSELoss()

for epoch in range(100):
    for batch in train_loader:
        optimizer.zero_grad()
        coarse_features = batch['coarse']  # Layer-avg quantities
        fine_truth = batch['fine']         # High-res Ri profile
        
        pred = model(coarse_features)
        loss = criterion(pred, fine_truth)
        loss.backward()
        optimizer.step()
```

---

### 2.2 Learning Effective φ Functions

**Goal:** Find grid-dependent $\phi_m^*(\zeta, \Delta z)$ and $\phi_h^*(\zeta, \Delta z)$ that reproduce bulk fluxes better than standard MOST.

**Method:** Symbolic regression (PySR) or neural ODE.

**Example (PySR — Automated Equation Discovery):**
```python
from pysr import PySRRegressor

# Training data: (zeta, delta_z) → phi_m_effective
X_train = df[['zeta_geom', 'delta_z_norm']]  # Normalized delta_z
y_train = df['phi_m_corrected']  # Effective value from validation

model_phi = PySRRegressor(
    niterations=50,
    binary_operators=["+", "*", "/"],
    unary_operators=["exp", "log"],
    constraints={'exp': 1, 'log': 1},  # Limit complexity
)
model_phi.fit(X_train, y_train)

# Best discovered equation
print(model_phi.equations_)
# Example output: phi_m_eff = 1 + 4.7*zeta + 0.3*log(delta_z/10)
```

---

## 3. Student Projects (Ordered by Difficulty)

### Project 1: Bias Ratio Predictor (Easiest)
- **Time:** 1-2 weeks
- **Goal:** Train Random Forest to predict $B$ from layer properties
- **Deliverable:** Jupyter notebook with MAE < 0.05

**Success Metric:** Model predicts $B$ within ±5% of analytic truth on test set.

### Project 2: Layer Classifier (Medium)
- **Time:** 2-3 weeks
- **Goal:** Classify layers as "needs correction" vs "skip"
- **Deliverable:** Confusion matrix showing >90% accuracy

**Success Metric:** False negatives (missed correction) < 5%.

### Project 3: Parameter Tuner (Advanced)
- **Time:** 3-4 weeks
- **Goal:** Learn $(\alpha, \gamma)$ from regime characteristics
- **Deliverable:** Meta-model + validation on 3 LES cases

**Success Metric:** Correction with ML-tuned params reduces flux RMSE by ≥20% vs fixed params.

### Project 4: Fine-Grid Emulator (Research-Level)
- **Time:** 4-6 weeks
- **Goal:** Neural net reconstructing high-res $Ri_g(z)$ from coarse inputs
- **Deliverable:** Conference-ready figure showing emulated vs LES profiles

**Success Metric:** Profile RMSE < 0.02 in stable regime.

---

## 4. Integration with McNider Correction Workflow

**Step-by-Step:**

```python
# 1. Load pre-trained ML models (done once at model initialization)
model_B = load_pickle('bias_predictor.pkl')
clf_layer = load_pickle('layer_classifier.pkl')

# 2. Runtime (every timestep, every layer)
for k in range(nz):
    # Extract layer properties
    delta_z = z[k+1] - z[k]
    z_geom = np.sqrt(z[k] * z[k+1])
    zeta = z_geom / L
    Ri_b_raw = compute_bulk_ri(k)  # From model state
    
    # ML Step 1: Should we correct this layer?
    needs_correction = should_apply_correction(
        clf_layer, B_pred=1.0, zeta=zeta, 
        Ri_b=Ri_b_raw, delta_z=delta_z
    )
    
    if not needs_correction:
        K_corrected[k] = K_original[k]
        continue
    
    # ML Step 2: Fast bias estimation
    B_ml = predict_bias_ratio(
        model_B, delta_z, z[k], z0, Ri_b_raw, 
        zeta, dtheta_dz[k], du_dz[k]
    )
    
    # ML Step 3: (Optional) Adaptive parameter tuning
    alpha_ml = predict_alpha(regime_features)  # If available
    
    # Physics Step: McNider correction
    f_c = compute_correction_factor(
        B_ml, zeta, delta_z, alpha=alpha_ml
    )
    K_corrected[k] = K_original[k] * f_c
```

**Computational Cost:**
- ML inference: ~0.1 ms per layer (Random Forest on CPU)
- Total overhead: <1% for typical 50-layer model
- **Trade-off:** Replaces expensive analytic integration (5-10 ms) with fast ML (0.1 ms) → net speedup!

---

## 5. Validation Checklist

Before deploying ML components:

- [ ] **Physical bounds:** All ML outputs respect $1.0 \leq B \leq 3.0$, $0.5 \leq \alpha \leq 2.5$
- [ ] **Neutral preservation:** ML-corrected profiles satisfy $Ri_g \to \zeta$ as $\zeta \to 0$
- [ ] **Generalization:** Test on unseen LES cases (GABLS2, GABLS3)
- [ ] **Robustness:** Inject 10% noise into features; check output stability
- [ ] **Interpretability:** Feature importance aligns with physics (e.g., $\Delta z$ and $\zeta$ dominate)

---

## 6. Code Repository Structure

```
ABL/
├── ML/
│   ├── models/
│   │   ├── bias_predictor.pkl       # Trained Random Forest for B
│   │   ├── layer_classifier.pkl     # Logistic Regression for flagging
│   │   ├── param_tuner.pkl          # Alpha/gamma meta-model
│   │   └── emulator_ckpt.pt         # Neural net for fine-grid reconstruction
│   ├── training/
│   │   ├── generate_synthetic_data.py
│   │   ├── train_bias_predictor.py
│   │   ├── train_classifier.py
│   │   └── train_emulator.py
│   ├── deployment/
│   │   ├── ml_inference.py          # Runtime ML functions
│   │   └── integration_example.py   # Full correction workflow
│   └── notebooks/
│       ├── 01_Synthetic_Data_Demo.ipynb
│       ├── 02_Bias_Predictor_Training.ipynb
│       ├── 03_Layer_Classification.ipynb
│       └── 04_Parameter_Tuning.ipynb
```

---

## 7. References & Further Reading

**ML for Atmospheric Physics:**
- Rasp et al. (2018): "Deep learning to represent subgrid processes in climate models"
- Reichstein et al. (2019): "Deep learning and process understanding for data-driven Earth system science"

**Richardson Number & SBL:**
- England & McNider (1995): Stability functions based upon shear functions
- Businger et al. (1971): Flux-profile relationships (MOST foundation)

**Practical ML:**
- Scikit-learn docs: https://scikit-learn.org
- PyTorch tutorials: https://pytorch.org/tutorials
- PySR (symbolic regression): https://github.com/MilesCranmer/PySR

---

## 8. Summary Table: ML Tools by Use Case

| Use Case | ML Method | Training Data | Output | Complexity |
|----------|-----------|---------------|--------|------------|
| **Fast B prediction** | Random Forest | Synthetic MOST profiles | Bias ratio B | ⭐ Easy |
| **Layer flagging** | Logistic Regression | Labeled synthetic profiles | Binary (correct/skip) | ⭐ Easy |
| **Parameter tuning** | RF Regressor | LES validation runs | (α, γ) optimal | ⭐⭐ Medium |
| **QC / Outlier detection** | Isolation Forest | Clean tower data | Anomaly flags | ⭐⭐ Medium |
| **Gradient imputation** | Gaussian Process | Tower time series | Missing ∂θ/∂z | ⭐⭐⭐ Advanced |
| **Fine-grid emulation** | Neural Network | Paired coarse/fine runs | High-res Ri profile | ⭐⭐⭐ Advanced |
| **φ function discovery** | Symbolic Regression | LES + tower ensemble | New φ equation | ⭐⭐⭐⭐ Research |

---

# ML Guide: Curvature-Aware SBL Corrections (practical)

Scope
- Train lightweight ML surrogates for:
  1) Ĝ(Δz, ζ, Ri_b, z_g, a_m,a_h) — multiplicative damping on K.
  2) P̂_laminar(Δz, ζ, Ri_b, Δz/z_g) — laminar/turbulent classifier.
  3) Ri_ĉ*(state) — dynamic critical Richardson (see Dynamic_Ric_ML_Strategy.md).

Feature engineering
- Core numeric: Δz, ζ, Ri_b, z_g, Δz/z_g, local gradients (dU/dz,dθ/dz), a_m,a_h proxies.
- Interaction features: ζ^2, (Δz/Δz_ref)^p, ζ*(Δz/Δz_ref).
- Physics-features: sign(H_sfc), u_*, θ_*, time-of-day.
- Normalize ζ by a reference ζ_ref for stability in training.

Model choices
- G surrogate: GradientBoostingRegressor / RandomForest with monotone post-processing and clipping to (0,1].
- Classifier: LogisticRegression or LightGBM with class weights; prefer probabilistic outputs.
- Ri_c*: Symbolic regression → simple closed form; keep coefficients small for stability.

Constraints & loss
- Enforce G(0)=1 and ∂G/∂ζ|_0≈0 via:
  - Augment training set with ζ≈0 examples and loss penalty on (G_pred−1)^2.
  - Use feature transform (zeta^2) so small-ζ slope ≈0.
- Clip G ∈ (ε,1]. Use ε ~ 1e-3 to avoid zeroing K in solver.
- Penalize negative monotonicity violations if physics requires monotone decrease.

Training recipe
1. Build synthetic dataset from analytic φ forms over grid sweep (Δz ∈ [5,200]m, ζ ∈ [0,1], a_m,a_h range).
2. Augment with LES/tower label matching when available.
3. Split by site / experiment to avoid leakage.
4. Train with early stopping and simple model selection (speed ≈ goal).
5. Export top model to ONNX + CSV LUT for fallback.

Evaluation (report)
- Regression: RMSE, MAE, bias near ζ→0 (target < small ε).
- Classification: AUC‑ROC, AUC‑PR, F1 at operational threshold (e.g., 0.999).
- Operational: change in bias ratio B before/after (target ≥ 40% reduction).
- Runtime: µs per call target (benchmark on intended CPU/Fortran host).

Deployment pattern
- Offer both: deterministic LUT (grid: Δz×ζ×Ri_b) and ONNX small model.
- Provide analytic fallback G_template for safety (see notebook).
- Version artifacts and log ModelVersion/DatasetVersion in diagnostics.

Notes
- Keep models tiny and interpretable for operational adoption.
- Document training data provenance and calibration hyperparameters.
