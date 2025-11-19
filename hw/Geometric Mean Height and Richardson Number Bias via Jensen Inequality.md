# Homework: Geometric Mean Height & Richardson Number Bias (Jensen's Inequality)

Course: Boundary-Layer Meteorology / Advanced Atmospheric Physics  
Level: Undergraduate and above — concise, runnable, and reproducible.

---

## Goal (one sentence)
Show that for concave-down gradient Richardson profiles typical of stable boundary layers, the layer-averaged bulk Richardson number underestimates the point (gradient) Richardson number at the geometric-mean height; quantify the bias and suggest practical corrections.

Learning objectives
- Apply Jensen's inequality to atmospheric profiles.
- Understand why the geometric mean z_g = √(z0 z1) is a good representative height for log-like profiles.
- Compute Ri_g and Ri_b from discrete data, compare, and propose a simple correction.

---

## Problem statement (short)
Given a layer [z0,z1], show that if Ri_g(z) is concave-down on that interval then
$$
Ri_b=\frac{1}{\Delta z}\int_{z_0}^{z_1} Ri_g(z)\,dz < Ri_g(z_g),\qquad z_g=\sqrt{z_0 z_1}.
$$
Provide a concise proof, a numerical verification, and a recommended pragmatic correction for models that cannot refine Δz.

---

## Part A — Math warm‑up (Jensen)
1. State Jensen's inequality for a concave function f on [a,b]:
   $$f\!\Big(\frac{1}{b-a}\int_a^b x\,dx\Big) \ge \frac{1}{b-a}\int_a^b f(x)\,dx.$$
2. Show $\ln z$ is concave for z>0: compute $\dfrac{d^2}{dz^2}\ln z=-1/z^2<0$.

Short answer expected: 2–3 lines.

---

## Part B — Why z_g (sketch)
1. For log wind $U(z)=(u_*/\kappa)\ln(z/z_0)+C$, show the exact layer-mean gradient equals the point gradient at the logarithmic mean
   $$z_{ \ln}=\frac{z_1-z_0}{\ln(z_1)-\ln(z_0)}.$$
2. Note ln z_g = (ln z0 + ln z1)/2, so z_g is the midpoint in log-space.
Short answer expected: 3–5 lines.

---

## Part C — Main proof (concise)
1. Apply Jensen to concave Ri_g(z): conclude
   $$Ri_b < Ri_g\!\left(\frac{z_0+z_1}{2}\right).$$
2. Use concavity of ln to show $z_g \le (z_0+z_1)/2$.
3. If Ri_g decreases with z (typical SBL), then
   $$Ri_g(z_g) \ge Ri_g\!\left(\frac{z_0+z_1}{2}\right) > Ri_b.$$
One paragraph, with the chain of inequalities.

---

## Part D — Numerical verification (required)
Analytic quadratic case: $Ri_g(z)=c_1 z + c_2 z^2$.

Run the snippet below, report z_g, Ri_g(z_g), Ri_b, and bias ratio B = Ri_g(z_g)/Ri_b. Include one plot (Ri_g vs z with shaded [z0,z1] and marker at z_g).

```python
# quick check: quadratic example
import math
z0, z1 = 10.0, 100.0
c1, c2 = 0.01, -0.0001
zg = math.sqrt(z0*z1)
Rig_zg = c1*zg + c2*zg*zg
Dz = z1 - z0
Rb = (0.5*c1*(z1*z1 - z0*z0) + (c2/3.0)*(z1**3 - z0**3))/Dz
B = Rig_zg / Rb
print(f"z_g={zg:.2f} m, Ri_g(z_g)={Rig_zg:.4f}, Ri_b={Rb:.4f}, B={B:.3f}")
```

Expected: Ri_b < Ri_g(z_g), B>1 (report values).

Deliverable: single Jupyter notebook cell with this code, a short figure, and one-paragraph interpretation.

---

## Part E — Practical correction (one-liner)
If you cannot refine Δz, evaluate stability at z_g instead of using Ri_b, or apply a multiplicative damping to diffusivity:
$$K^* = K\cdot \exp\!\Big[-D\big(\tfrac{\Delta z}{\Delta z_r}\big)^p\big(\tfrac{\zeta}{\zeta_r}\big)^q\Big],\quad q\ge2$$
State chosen constants and justify briefly.

---

## Deliverables (clear)
1. PDF or Markdown report (≤3 pages) with proofs (Parts A–C), numeric results (Part D), figure, and correction recommendation (Part E).
2. Jupyter notebook (runnable) with the snippet and plot.
3. One-slide summary (1 PNG) of results and practical recommendation.

---

## Grading rubric (100 pts)
- Correct math / proof: 25
- Numerical implementation & correctness: 30
- Figures and interpretation (clear axis/legend): 20
- Reproducibility (notebook runs): 15
- Concise policy/practical recommendation: 10

---

## Hints & Resources (quick)
- Use z_g = sqrt(z0*z1) for representative height in log-like profiles.
- For gradients from discrete data, use centered differences for interior points; forward/backward at boundaries.
- Prefer Simpson/trapezoid integration for Ri_b if you have 3+ points in the layer.
- Suggested datasets: ARM SGP, CASES-99, Cabauw, GABLS LES.
- Python libs: numpy, matplotlib, xarray (if using NetCDF).

---

## Extension (optional, +15%)
1. Compute Simpson and trapezoid Ri_b and report their relative errors vs analytic Ri_b for the quadratic example.
2. Implement the simple K damping above and show how K at z_g changes (percent) for D tuned to reduce B by ~40% at Δz=60 m.

---

## Submission
- Upload notebook (.ipynb) and report (.pdf or .md) to course repo or LMS.
- Filename convention: `hw01_lastname.ipynb` and `hw01_lastname.pdf`.
- Due: [instructor sets date].

Good luck — keep answers concise and reproducible.

