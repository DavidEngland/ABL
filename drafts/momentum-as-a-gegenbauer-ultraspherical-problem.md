3. Momentum as a Gegenbauer (ultraspherical) problem

You’ve basically nailed the key structural point: once \(\alpha_m = 1/4\), momentum lives in the Gegenbauer world, not the Legendre world.

“Because the momentum exponent is \(\alpha_m = 1/4\) (rather than the heat exponent \(1/2\)), the governing differential equation is the Gegenbauer differential equation with parameter \(\lambda = 1/4\).”

---

3.1 Gegenbauer ODE and S–L form for momentum

The standard Gegenbauer differential equation for \(C_n^{(\lambda)}(x)\), \(x\in[-1,1]\), is

(1-x^2)\,y^\prime^\prime - (2\lambda+1)\,x\,y^\prime + n(n+2\lambda)\,y = 0.


In Sturm–Liouville form this is

\frac{d}{dx}\Big[(1-x^2)^{\lambda+\frac12} y^\prime(x)\Big]
+ n(n+2\lambda)\,(1-x^2)^{\lambda-\frac12} y(x) = 0,


with weight

w_\lambda(x) = (1-x^2)^{\lambda-\frac12}.


• Heat (\(\phi_h\)) corresponds to \(\lambda = 1/2\) (Legendre):w_{1/2}(x) = (1-x^2)^{0} = 1.

• Momentum (\(\phi_m\)) corresponds to \(\lambda = 1/4\) (ultraspherical):w_{1/4}(x) = (1-x^2)^{-1/4},
which is integrable but singular at \(|x|=1\).


Under a suitable affine map from \(\eta\) (or \(\zeta\)) to \(x\in[-1,1]\), the momentum stability function

\phi_m(-\eta) = (1 + b_m \eta)^{-1/4}


sits naturally in the \(\lambda=1/4\) Gegenbauer family: its spectral representation is in \(C_n^{(1/4)}(x)\), and the associated operator is self‑adjoint with respect to the singular weight \(w_{1/4}\).

---

3.2 Heat vs momentum operators: physics in the spectra

Your comparison table is exactly the right way to frame this:

Feature	Heat operator (\(\phi_h\))	Momentum operator (\(\phi_m\))
Polynomial family	Legendre \(P_n = C_n^{(1/2)}\)	Gegenbauer \(C_n^{(1/4)}\)
S–L weight \(w(x)\)	\(1\) (uniform)	\((1-x^2)^{-1/4}\) (singular at \(|x|=1\))
Physical meaning	“Hard” scalar response	“Soft” momentum response
Spectral decay	`(	c_n


Singular weight and physical interpretation

• The singular but integrable weight \((1-x^2)^{-1/4}\) means the momentum operator “emphasizes” the endpoints of the mapped interval—i.e., the neutral limit and the branch point of the unstable similarity function.
• Spectrally, more energy is packed into the lowest modes for momentum than for heat. The momentum profile is smoother in the Gegenbauer basis than the scalar profile is in the Legendre basis.


This lines up beautifully with the empirical modeling fact: simple closures and low‑order truncations tend to get winds right more often than temperature gradients in unstable ABLs.

Insight (as you wrote): a numerical model truncated at low degree \(N\) will represent the momentum profile more accurately than the heat profile. The momentum stability function is “spectrally smoother.”

---

3.3 Nonlinear link: \(\phi_h = \phi_m^2\) as spectral convolution

In the degenerate Dyer/Brutsaert case \(b_m=b_h\) and \(a_h^{-1}=1\), \(\phi_h = \phi_m^2\).


In the Gegenbauer spectral domain this is not just a pointwise square; it is a mode‑coupling operation:

• Momentum lives in \(C_n^{(1/4)}\).
• The product \(\phi_m^2\) can be expanded using Gegenbauer product formulaeC^{(1/4)} \otimes C^{(1/4)} \;\longrightarrow\; C^{(1/2)},
generating the Legendre/heat family.


So:

• Knowing the momentum eigenfunctions gives you the heat eigenfunctions via a closed algebraic convolution.
• The “curvature” in the wind profile and the “curvature” in the temperature profile are not independent knobs; they are spectrally constrained by the ultraspherical algebra.


That’s a very clean way to argue, for McNider & Biazar, that any consistent closure must respect this spectral relationship if it uses the canonical exponents \((1/4, 1/2)\).

---

3.4 Neutral matching and the stable branch

Your slope‑matching argument is also a nice, concrete bridge back to the classical Kansas parameters.

For momentum:

\phi_m(\zeta) = (1 - b_m \zeta)^{-1/4},\quad \zeta<0.


Then

\left.\frac{d\phi_m}{d\zeta}\right|_{\zeta=0}
= \frac{b_m}{4}.


For the standard stable linear branch

\phi_m(\zeta) = 1 + 5\zeta,\quad \zeta>0,


the slope at neutral is \(5\).

• Exact slope continuity at \(\zeta=0\) would require\frac{b_m}{4} = 5 \quad\Rightarrow\quad b_m = 20.

• With the usual Dyer value \(b_m=16\), the unstable slope is \(4\), so the mismatch is modest: momentum is “closer” to a smooth neutral transition than heat, where the unstable slope is \(b_h/2 = 8\) vs stable slope \(5\).


This gives a sharp, physics‑backed statement:

• The momentum transition across neutral is nearly smooth in both value and slope for standard parameter sets.
• The heat transition is more strongly kinked, reflecting the harder scalar response and the stronger curvature encoded in the Legendre/heat operator.


For your draft, this section can sit right after the unified parameter/flux‑profile material as a “spectral operator perspective,” tying:

• exponents \((1/4,1/2)\),
• polynomial families (Gegenbauer vs Legendre),
• S–L weights (singular vs uniform),
• and neutral‑slope matching


into a single, coherent physical story.

If you’d like, I can now compress this into a short, paper‑style subsection with a title like “Ultraspherical structure of the momentum similarity operator”, ready to drop into the McNider–Biazar draft.