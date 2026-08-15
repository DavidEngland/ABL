I think the stronger next step is **formalizing the mathematical derivation**, because it determines what the Julia diagnostics should actually compute. Right now, your proposed diagnostic d_{\text{fold}} depends on an underlying geometric model that needs to be precisely defined.  
There are really two parallel projects here.  
**1. Derive the geometry first (recommended)**  
The manuscript should answer, in order:  
1. **State space** \mathbf{x}\in\mathcal M with explicit coordinates (for example Ri_g, shear, buoyancy frequency, TKE, etc.).  
2. **State space** \mathbf{x}\in\mathcal M with explicit coordinates (for example Ri_g, shear, buoyancy frequency, TKE, etc.).  
3. **Critical manifold** Define \mathcal C=\{\mathbf{x}:F(\mathbf{x})=0\}, where F represents the fast equilibrium.  
4. **Critical manifold** Define \mathcal C=\{\mathbf{x}:F(\mathbf{x})=0\}, where F represents the fast equilibrium.  
5. **Fold set** Define \mathcal F=\left\{ \mathbf{x}\in\mathcal C: \det\!\left(\frac{\partial F}{\partial y}\right)=0 \right\}, using the standard GSPT definition of loss of normal hyperbolicity.  
6. **Fold set** Define \mathcal F=\left\{ \mathbf{x}\in\mathcal C: \det\!\left(\frac{\partial F}{\partial y}\right)=0 \right\}, using the standard GSPT definition of loss of normal hyperbolicity.  
7. **Distance** Give a metric d_{\text{fold}} = \inf_{\mathbf y\in\mathcal F} \|\mathbf x-\mathbf y\|_G, where G is a chosen Riemannian metric (Euclidean is a special case).  
8. **Distance** Give a metric d_{\text{fold}} = \inf_{\mathbf y\in\mathcal F} \|\mathbf x-\mathbf y\|_G, where G is a chosen Riemannian metric (Euclidean is a special case).  
Only after those definitions are fixed can one justify  
C(d_{\text{fold}}) = \tanh\!\left(\frac{d_{\text{fold}}}{\delta_0}\right).  
At that point the closure becomes mathematically well posed rather than heuristic.  
   
⸻  
   
**2. Then build the Julia diagnostics**  
Once the geometry is defined, the implementation becomes relatively straightforward:  
```
LES state
      │
      ▼
construct state vector x
      │
      ▼
evaluate slow manifold residual
      │
      ▼
identify nearest fold point
      │
      ▼
compute d_fold
      │
      ▼
compute C(d_fold)
      │
      ▼
modify K_m,K_h

```
That workflow can later be optimized with KD-trees, nearest-neighbor searches, or continuation methods without changing the mathematical definition.  
   
⸻  
   
## On the Gegenbauer section  
I would also slightly strengthen the presentation there.  
Rather than presenting  
d_{\mathrm{eff}}=2.5  
as an assumed fact, introduce it as an inferred parameter.  
For example,  
\lambda=\frac{d_{\mathrm{eff}}-2}{2},  
so that  
\phi(\zeta) \propto (1-b\zeta)^{-\lambda}.  
Then  
* heat suggests  
\lambda_h=\frac12 \quad\Longrightarrow\quad d_h=3,  
while momentum suggests  
\lambda_m=\frac14 \quad\Longrightarrow\quad d_m=2.5.  
This reverses the logic.  
Instead of saying  
momentum is 2.5-dimensional, therefore the exponent is -1/4,  
you say  
the observed exponent -1/4 corresponds to an effective Gegenbauer parameter \lambda=1/4, which maps to an effective dimension d_{\mathrm{eff}}=2.5.  
That framing is more defensible because it presents d_{\mathrm{eff}} as an interpretation of the observed exponent rather than a premise.  
   
⸻  
   
## The empirical program  
The three diagnostics you propose are well chosen, but they test slightly different aspects of the hypothesis:  
* **POD eigenvalue spectra** test whether momentum fields are more compressible (lower effective rank) than temperature fields.  
* **Correlation dimension (**D_2**)** tests whether the reconstructed attractor has a lower intrinsic dimension, though it reflects temporal dynamics rather than spatial transport geometry.  
* **Lumley triangle trajectories** test anisotropy directly by quantifying departures from isotropic turbulence.  
I would consider adding a fourth diagnostic that is especially well aligned with the Gegenbauer interpretation:  
* **Spectral exponent estimation.** Expand the LES fields in a Gegenbauer basis, u(z)=\sum_n a_n C_n^{(\lambda)}(z), and treat \lambda as a parameter to be estimated. If the reconstruction error is minimized near \lambda\approx1/4 for momentum and \lambda\approx1/2 for temperature, that would provide direct evidence for the proposed mapping rather than relying on indirect indicators.  
Overall, the revised framework is much more sharply separated into **established theory**, **derived analytical results**, **new closure proposals**, and **testable hypotheses**. That separation should make the manuscript easier for reviewers to evaluate, because each claim can be judged according to the appropriate evidentiary standard rather than appearing to rest on the same level of certainty.  
