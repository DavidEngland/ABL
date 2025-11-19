# Homework: Geometric Mean Height & Richardson Number Bias (Jensen's Inequality)

Course: Boundary-Layer Meteorology / Advanced Atmospheric Physics  
Level: Undergraduate and above

---

## Problem statement (brief)
Atmospheric models compute a bulk Richardson number across a layer,
$$
Ri_b = \frac{g}{\theta}\frac{\Delta\theta\,\Delta z}{(\Delta U)^2},
$$
while the local (gradient) Richardson number is
$$
Ri_g(z)=\frac{(g/\theta)\,\partial\theta/\partial z}{(\partial U/\partial z)^2}.
$$
Show that when \(Ri_g(z)\) is concave-down over \([z_0,z_1]\) (typical SBL), the layer average underestimates the point value at the geometric-mean height \(z_g=\sqrt{z_0 z_1}\):
$$
Ri_b = \frac{1}{\Delta z}\int_{z_0}^{z_1} Ri_g(z)\,dz < Ri_g(z_g).
$$

---

## Part A — Jensen's inequality & concavity (warm-up)

A1. State Jensen's inequality for a concave function \(f\):
For concave \(f\) on \([a,b]\),
$$
f\!\left(\frac{1}{b-a}\int_a^b x\,dx\right) \ge \frac{1}{b-a}\int_a^b f(x)\,dx.
$$
Concavity (twice differentiable): \(f''(x)\le 0\) on \([a,b]\).

A2. Show \(\ln z\) is concave for \(z>0\):
$$
\frac{d^2}{dz^2}\ln z = -\frac{1}{z^2}<0,\quad z>0.
$$

---

## Part B — Why the geometric mean?

B1. Logarithmic mean for log wind:
For \(U(z)=\frac{u_*}{\kappa}\ln(z/z_0)+C\),
the exact layer-averaged gradient equals the point gradient at the logarithmic mean
$$
\bar z_{\ln}=\frac{z_1-z_0}{\ln(z_1)-\ln(z_0)}.
$$

B2. Geometric mean is midpoint in log space:
$$
\ln z_g=\tfrac12(\ln z_0+\ln z_1)\quad\Rightarrow\quad z_g=\sqrt{z_0 z_1}.
$$

B3. For thin layers (\(z_1/z_0\to1\)) the log-mean and geometric mean coincide to \(O((\ln(z_1/z_0))^2)\); hence \(z_g\) is a practical representative height for log-like profiles.

---

## Part C — Concave-down \(Ri_g\) and bias proof

C1. Apply Jensen to concave \(Ri_g(z)\):
Since \(Ri_g\) concave → by Jensen
$$
\frac{1}{\Delta z}\int_{z_0}^{z_1} Ri_g(z)\,dz \le Ri_g\!\left(\frac{1}{\Delta z}\int_{z_0}^{z_1} z\,dz\right) = Ri_g\!\left(\frac{z_0+z_1}{2}\right).
$$
So \(Ri_b < Ri_g\) at the arithmetic mean.

C2. Use concavity of \(\ln z\) to compare means:
Because \(\ln\) is concave,
$$
\ln z_g = \tfrac12(\ln z_0+\ln z_1) > \ln\!\left(\tfrac{z_0+z_1}{2}\right),
$$
hence \(z_g < (z_0+z_1)/2\). For typical SBL where \(Ri_g\) decreases with height (near-surface), this implies
$$
Ri_g(z_g) > Ri_g\!\left(\tfrac{z_0+z_1}{2}\right),
$$
so \(Ri_b < Ri_g(z_g)\).

C3. Alternative, cleaner route: change variable \(s=\ln z\). If \(\tilde Ri(s)=Ri_g(e^s)\) is concave in \(s\), apply Jensen on \(s\)-interval to get
$$
\frac{1}{\ln(z_1/z_0)}\int_{\ln z_0}^{\ln z_1}\tilde Ri(s)\,ds < \tilde Ri\!\left(\tfrac{\ln z_0+\ln z_1}{2}\right)=Ri_g(z_g).
$$
For thin layers the log-weighted average approximates the arithmetic layer average, yielding \(Ri_b \lesssim Ri_g(z_g)\).

---

## Part D — Numerical verification (analytic quadratic example)

Assume \(Ri_g(z)=c_1 z + c_2 z^2\) with \(c_1>0, c_2<0\).
Example: \(z_0=10\) m, \(z_1=100\) m, \(c_1=0.01\), \(c_2=-0.0001\).
Compute:
- \(z_g=\sqrt{z_0 z_1}\).
- \(Ri_g(z_g)=c_1 z_g + c_2 z_g^2\).
- Analytic integral:
$$
Ri_b=\frac{1}{\Delta z}\int_{z_0}^{z_1}(c_1 z + c_2 z^2)\,dz
=\frac{1}{\Delta z}\left[\tfrac{c_1}{2}(z_1^2-z_0^2)+\tfrac{c_2}{3}(z_1^3-z_0^3)\right].
$$

Short Python snippet to reproduce:
```python
# compute example values
import math
z0, z1 = 10.0, 100.0
c1, c2 = 0.01, -0.0001
zg = math.sqrt(z0*z1)
Rig_zg = c1*zg + c2*zg*zg
Dz = z1 - z0
Rb = (0.5*c1*(z1*z1 - z0*z0) + (c2/3.0)*(z1**3 - z0**3))/Dz
print("z_g =", zg)
print("Ri_g(z_g) =", Rig_zg)
print("Ri_b =", Rb)
```
You should observe `Ri_b < Ri_g(z_g)` for the example.

---

## Part E — Interpretation and practical correction

E1. Why underestimation leads to overmixing: models map \(Ri\) → eddy diffusivities \(K_m,K_h\) (decreasing functions of \(Ri\)). If \(Ri_b\) underestimates true local stability, computed \(K\) are too large → excessive vertical mixing.

E2. Using arithmetic mean height worsens the bias (since \(z_g<z_{\text{arith}}\) and \(Ri_g\) typically decreases with height in SBL).

E3. Practical correction (no grid change): evaluate representative Ri at \(z_g\) (or use log-mean) when computing stability functions; or apply a curvature-aware multiplicative correction \(G(\zeta,\Delta z)\) to \(K\) (ensure \(G(0)=1\), \(G'(0)=0\)).

---

## Optional challenge (derivation hint)

Taylor expand \(Ri_g\) about \(z_g\):
$$
Ri_g(z)\approx Ri_g(z_g) + Ri_g'(z_g)(z-z_g) + \tfrac12 Ri_g''(z_g)(z-z_g)^2.
$$
Integrate over \([z_0,z_1]\); the linear term integrates to zero about the symmetric point \(z_g\), leaving a quadratic bias term proportional to \(Ri_g''(z_g)\). Use this to construct a leading-order correction factor \(G\).

---

## Submission
- Format: PDF or Markdown with embedded math (LaTeX rendered).
- Include code used for numerical checks.
- Provide figures where helpful (plot \(Ri_g\), shaded layer, indicate \(z_g\), show area under curve vs point value).

Good luck.
