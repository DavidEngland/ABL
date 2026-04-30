# spectral_scalar_closure.py

## Purpose

This module turns the recent Gegenbauer / Legendre analysis into an executable closure model for momentum, heat, and other scalars in the atmospheric boundary layer (ABL).

The main idea is simple:

1. Represent the **momentum similarity function** in a bounded spectral coordinate.
2. Obtain **heat or other scalar similarity functions** by applying a mode-dependent filter to the same momentum coefficients.
3. Derive **Prandtl number**, **Richardson number**, and **other scalar-transfer ratios** from that shared representation.

So the module is not doing an arbitrary fit for each variable. It is enforcing a common dynamical backbone and then letting the scalar transport differ through an interpretable filter scale.

---

## Mathematical Model

### 1. Compactify stability space

Instead of expanding directly in $\zeta = z/L$, the code maps stability to a bounded coordinate:

$$
\xi(\zeta) = \tanh(\alpha \zeta), \qquad \xi \in [-1,1].
$$

Why this matters:

- $\zeta$ is unbounded on the stable side.
- Gegenbauer and Legendre polynomials naturally live on $[-1,1]$.
- The parameter $\alpha$ controls how quickly the map moves from near-neutral to strongly stable conditions.

### 2. Expand momentum in Gegenbauer modes

The momentum similarity function is represented as

$$
\phi_m(\zeta) = \sum_{n=0}^{N} a_n \, C_n^{(\lambda)}\!\bigl(\xi(\zeta)\bigr),
$$

where:

- $a_n$ are the momentum coefficients,
- $C_n^{(\lambda)}$ are Gegenbauer polynomials,
- $\lambda$ controls the basis weighting.

Special cases:

- $\lambda = 1/2$ gives the **Legendre** basis.
- $\lambda = 1/4$ is the canonical unstable momentum generating-function family.

### 3. Generate heat or other scalars by filtering the same modes

For a scalar $s$ (heat, humidity, water vapor, $CO_2$, etc.), the code uses

$$
\phi_s(\zeta) = S_{0,s} \sum_{n=0}^{N} r_{n,s} \, a_n \, C_n^{(\lambda)}\!\bigl(\xi(\zeta)\bigr),
$$

with

$$
r_{n,s} = \exp\!\left(-\frac{n}{n_{c,s}}\right).
$$

Interpretation:

- the scalar keeps the same large-scale momentum structure,
- but higher-order modes are damped more strongly,
- and the damping strength is controlled by $n_{c,s}$.

Small $n_c$ means aggressive filtering. Large $n_c$ means the scalar stays closer to momentum.

---

## Physical Interpretation

### Why this is better than separate curve fits

The code is enforcing this principle:

> momentum, heat, and other scalars are not independent closures; they are related transport channels built on the same structure.

That gives three advantages:

1. **Shared physics**: one momentum expansion supports multiple scalars.
2. **Interpretability**: scalar differences are encoded in $S_0$ and $n_c$, not hidden in unrelated coefficients.
3. **Consistency**: $Pr_t$, $Ri_g$, and scalar-transfer ratios are derived from the same objects.

### Hockey-stick curves

The “hockey-stick” behavior in $Pr_t(\zeta)$ comes from the fact that high-order heat modes are damped more strongly than momentum modes. Near neutral, the ratio is flat. In the transition range, the ratio bends upward. In very stable conditions, it approaches a plateau or slow-growth regime depending on the chosen coefficients and map.

---

## Connection to the Gegenbauer / Spherical Analysis

The code directly reflects the mathematical results developed in the notes.

### 1. Unstable branch generating functions

Canonical unstable MOST forms are generating functions of ultraspherical families:

$$
\phi_m \sim (1 - b_m \zeta)^{-1/4},
\qquad
\phi_h \sim (1 - b_h \zeta)^{-1/2}.
$$

These correspond to:

- Gegenbauer with $\lambda = 1/4$ for momentum,
- Legendre with $\lambda = 1/2$ for heat.

### 2. Exact anchor case

In the matched unstable case

$$
a_h = 1, \qquad b_m = b_h = 16,
$$

we have the exact identities

$$
\phi_h = \phi_m^2,
\qquad
Ri_g = \zeta,
\qquad
Pr_t = \phi_m.
$$

This is why the module includes `exact_ubl_anchor()`: it provides a benchmark solution that any generalized closure should recover in that limit.

### 3. What is known vs what must be fitted

Usually constrained by analysis:

- the shared momentum-basis structure,
- the use of Gegenbauer/Legendre families,
- the matched unstable anchor identities,
- the compactifying map shape.

Usually fitted from data:

- the neutral scalar ratio $S_0$,
- the filter scale $n_c$,
- site/regime dependence under very stable or intermittent conditions.

---

## Code Structure

## `SpectralScalarClosure`

This is the main class.

### Constructor

```python
SpectralScalarClosure(coeffs, lambda_=0.5, alpha=0.5)
```

Arguments:

- `coeffs`: the momentum coefficients $a_n$
- `lambda_`: Gegenbauer parameter
- `alpha`: width of the map $\xi(\zeta)=\tanh(\alpha\zeta)$

Validation:

- `coeffs` must be a non-empty 1D array
- `lambda_ > 0`
- `alpha > 0`

### `xi(zeta)`

Returns the bounded stability coordinate:

```python
xi = np.tanh(alpha * zeta)
```

### `_basis_matrix(zeta)`

Evaluates the Gegenbauer basis functions at all requested $\zeta$ values using SciPy’s `eval_gegenbauer`.

Output shape is essentially:

- rows = mode number $n$
- columns = evaluation points

### `phi_m(zeta)`

Builds momentum from the basis matrix and coefficient vector:

$$
\phi_m = \sum a_n C_n^{(\lambda)}.
$$

### `filter_weights(n_c)`

Returns the exponential modal weights:

$$
r_n = e^{-n/n_c}.
$$

This is the key “spectral damping” control.

### `phi_scalar(zeta, s0, n_c)`

Builds a scalar similarity function from the momentum coefficients plus filter weights:

$$
\phi_s = S_0 \sum r_n a_n C_n^{(\lambda)}.
$$

### `scalar_ratio(zeta, s0, n_c)`

Returns

$$
\frac{\phi_s}{\phi_m},
$$

which is the general scalar analog of turbulent Prandtl or Schmidt number.

### `prandtl(zeta, pr0=0.85, n_c=1.2)`

Convenience wrapper for heat:

$$
Pr_t = \frac{\phi_h}{\phi_m}.
$$

### `ri_g(zeta, pr0=0.85, n_c=1.2)`

Computes the gradient Richardson number from the internally generated $\phi_m$ and $\phi_h$:

$$
Ri_g = \zeta \frac{\phi_h}{\phi_m^2}.
$$

So once the model is specified, Richardson number is derived rather than independently fitted.

### `fit_scalar_filter(...)`

Fits either:

- both `s0` and `n_c`, or
- only `n_c` while holding `s0` fixed.

It uses `scipy.optimize.least_squares` to minimize the mismatch between observed scalar ratios and the model-generated ratio.

This is the main calibration method.

---

## Helper Functions

## `exact_ubl_anchor(zeta, b=16.0)`

Returns the exact matched unstable anchor solution:

- `phi_m = (1 - b zeta)^(-1/4)`
- `phi_h = phi_m**2`
- `pr_t = phi_m`
- `ri_g = zeta`

Use this as:

- a unit test,
- a benchmark,
- a calibration constraint.

## `recommend_calibration_strategy()`

Returns a dictionary describing:

- what is known from analysis,
- what needs to be fitted from data,
- where machine learning can help.

This is intended for notebook or workflow guidance, not physics itself.

---

## Typical Workflow

### 1. Choose a momentum basis

Start with a small coefficient vector, for example:

```python
coeffs = [1.0, 0.6, 0.2]
model = SpectralScalarClosure(coeffs=coeffs, lambda_=0.5, alpha=0.5)
```

### 2. Generate heat or scalar ratios

```python
zeta = np.array([0.0, 0.2, 0.5, 1.0, 2.0])
pr_t = model.prandtl(zeta, pr0=0.85, n_c=1.4)
```

### 3. Fit unknown constants from data

If you already know the neutral limit and only need the filter scale:

```python
fit = model.fit_scalar_filter(zeta, pr_t_obs, s0_init=0.85, n_c_init=1.0, fit_s0=False)
```

If both are unknown:

```python
fit = model.fit_scalar_filter(zeta, pr_t_obs, s0_init=0.85, n_c_init=1.0, fit_s0=True)
```

### 4. Extend to humidity or CO2

Use the same `coeffs`, `lambda_`, and `alpha`, but fit scalar-specific `s0` and `n_c`:

```python
sc_q = model.scalar_ratio(zeta, s0=s0_q, n_c=n_c_q)
sc_co2 = model.scalar_ratio(zeta, s0=s0_co2, n_c=n_c_co2)
```

### 5. Example: Methane tracer (CH4)

This is the exact pattern you outlined, using the tracer wrapper names:

```python
from spectral_scalar_closure import SpectralScalarClosure

# Initialize with momentum coefficients
model = SpectralScalarClosure(coeffs=my_momentum_coeffs, lambda_=0.25)

# Fit n_c and s0 from methane observations (ratio = phi_ch4 / phi_m)
fit_results = model.fit_tracer_filter(zeta_obs, methane_ratio_obs)

print(f"Neutral ratio: {fit_results['s0']}")
print(f"Filter scale (n_c): {fit_results['n_c']}")

# Predict methane transfer ratio on any zeta grid
methane_ratio_model = model.tracer_ratio(
	zeta_grid,
	s0=fit_results["s0"],
	n_c=fit_results["n_c"],
)
```

Notes:

- `methane_ratio_obs` should represent the transfer ratio $\phi_{CH4}/\phi_m$.
- If your neutral methane ratio is known from independent analysis, set `fit_s0=False` and pass that value as `s0_init`.
- Typical first pass is `lambda_=0.25` (momentum-like ultraspherical weighting) or `lambda_=0.5` (Legendre baseline), then compare fit residuals.

---

## What This Module Does Not Yet Do

This is a prototype, not a full production closure. It does **not** yet include:

- automatic estimation of the momentum coefficients from field data,
- explicit monotonicity or no-pole constrained Pad\'e surrogates,
- regime classification,
- WRF/Fortran integration,
- uncertainty quantification around fitted parameters.

Those are natural next steps.

---

## Should Unknown Constants Be Determined from Theory, Data, or ML?

Best answer: **all three, but in the right order**.

### Theory

Use theory to fix the architecture:

- basis family,
- map,
- anchor identities,
- shared momentum backbone.

### Data

Use real-world data to determine:

- neutral scalar ratios,
- filter scales,
- regime dependence.

### ML

Use ML only in a supporting role:

- predict $n_c$ from metadata,
- learn residual corrections,
- classify regimes.

Do not use ML as a fully unconstrained replacement for the spectral model.

---

## Recommended Next Step

The most useful immediate use of this module is:

1. fit heat data first to estimate $Pr_0$ and $n_c$,
2. test whether humidity and other scalars can share the same momentum basis,
3. use the exact UBL anchor case as a benchmark,
4. only then reduce the model to a Pad\'e or Fortran-ready implementation.

That keeps the mathematical structure intact while using observations only where genuine uncertainty remains.