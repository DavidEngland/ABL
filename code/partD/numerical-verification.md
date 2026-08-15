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