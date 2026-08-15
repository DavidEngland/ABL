Let’s work through this special case carefully and in full detail.

We take the model form for the gradient Richardson number:

Ri_g(\zeta) = \frac{\zeta}{1 + \beta \zeta}.


This is a rational function often used to mimic MOST‑like stability functions.

---

1. Series Expansion (small `\(\zeta\)`)

We expand in powers of \(\zeta\) using the geometric series:

Ri_g(\zeta) = \zeta \cdot \frac{1}{1+\beta \zeta}
= \zeta \left(1 - \beta \zeta + \beta^2 \zeta^2 - \beta^3 \zeta^3 + \beta^4 \zeta^4 - \cdots \right).


So:

Ri_g(\zeta) \approx \zeta - \beta \zeta^2 + \beta^2 \zeta^3 - \beta^3 \zeta^4 + \beta^4 \zeta^5 + \cdots


That’s the first 4–5 terms.

---

2. First and Second Derivatives

Ri_g^\prime(\zeta) = \frac{d}{d\zeta}\left(\frac{\zeta}{1+\beta \zeta}\right)
= \frac{(1+\beta \zeta) - \beta \zeta}{(1+\beta \zeta)^2}
= \frac{1}{(1+\beta \zeta)^2}.


Ri_g^\prime^\prime(\zeta) = \frac{d}{d\zeta}\left(\frac{1}{(1+\beta \zeta)^2}\right)
= -\frac{2\beta}{(1+\beta \zeta)^3}.


---

3. Curvature Interpretation

• \(Ri_g^\prime(\zeta)\) is always positive → the function increases monotonically.
• \(Ri_g^\prime^\prime(\zeta)\) is always negative (for \(\beta>0\)) → the function is concave‑down.
• This concavity is exactly what produces bulk vs gradient bias: the average slope (bulk Richardson) is smaller than the local slope, so \(B>1\).


---

4. Forming the Bias Ratio `\(B\)`

Suppose the bulk Richardson number is defined over a layer \([0,\zeta]\):

Ri_b(\zeta) = \frac{\Delta \theta \cdot \Delta z}{(\Delta U)^2} \;\;\; \text{(schematically)}.


In this toy model, we can approximate \(Ri_b\) as the layer average of \(Ri_g\):

Ri_b(\zeta) \approx \frac{1}{\zeta}\int_0^\zeta Ri_g(\xi)\,d\xi.


Compute:

Ri_b(\zeta) = \frac{1}{\zeta}\int_0^\zeta \frac{\xi}{1+\beta \xi}\,d\xi.


Integral:

\int \frac{\xi}{1+\beta \xi}\,d\xi = \frac{1}{\beta^2}\left(\beta \xi - \ln(1+\beta \xi)\right).


So:

Ri_b(\zeta) = \frac{1}{\zeta \beta^2}\left(\beta \zeta - \ln(1+\beta \zeta)\right).


Thus:

B(\zeta) = \frac{Ri_g(\zeta)}{Ri_b(\zeta)}
= \frac{\zeta/(1+\beta \zeta)}{\frac{1}{\zeta \beta^2}(\beta \zeta - \ln(1+\beta \zeta))}.


Simplify:

B(\zeta) = \frac{\zeta^2 \beta^2}{(1+\beta \zeta)(\beta \zeta - \ln(1+\beta \zeta))}.


This is the exact bias ratio for this special case.

---

5. Correction Factor `\(f_c\)`

ODE driver:

\frac{d \ln f_c}{d \ln \Delta z} = -\alpha (B-1)\left(\frac{\zeta}{\zeta_{\text{ref}}}\right)^q.


Solution:

f_c(\Delta z,\zeta) = \left(\frac{\Delta z}{\Delta z_{\text{ref}}}\right)^{-\alpha (B(\zeta)-1)(\zeta/\zeta_{\text{ref}})^q}.


So in this special case, the correction factor is explicitly tied to the rational form of \(Ri_g\) and the integral form of \(Ri_b\).

---

6. Dynamic Critical Richardson Number

Define:

Ri_{cr,\text{dyn}} = Ri_{cr,0}[1 + \gamma (B(\zeta)-1)].


• If curvature is strong (\(B\gg 1\)), the effective critical Richardson number is raised, allowing turbulence to persist longer.
• If curvature is weak (\(B\approx 1\)), the dynamic critical Richardson number reduces to the canonical value.


---

✅ Summary for Students

• Expansion: \(Ri_g = \zeta - \beta \zeta^2 + \beta^2 \zeta^3 - \cdots\).
• Derivatives: \(Ri_g^\prime = (1+\beta \zeta)^{-2}, \; Ri_g^\prime^\prime = -2\beta (1+\beta \zeta)^{-3}\).
• Curvature: concave‑down, so bulk underestimates stability.
• Bias ratio:B(\zeta) = \frac{\zeta^2 \beta^2}{(1+\beta \zeta)(\beta \zeta - \ln(1+\beta \zeta))}.

• Correction factor:f_c = \left(\frac{\Delta z}{\Delta z_{\text{ref}}}\right)^{-\alpha (B-1)(\zeta/\zeta_{\text{ref}})^q}.

• Dynamic critical Richardson number:Ri_{cr,\text{dyn}} = Ri_{cr,0}[1+\gamma(B-1)].



---

Would you like me to numerically evaluate \(B(\zeta)\) for a few values of \(\beta\) and \(\zeta\) (say \(\beta=1, \zeta=0.5,1.0,2.0\)) so you can see how curvature bias grows and how it feeds into \(f_c\) and \(Ri_{cr,\text{dyn}}\)?