For the GeoABL/GSPT framework you’ve been developing, these approaches are not equally well aligned with the underlying physics.  
**The strongest hierarchy is:**  
**The strongest hierarchy is:**  
1. **Replace the stability functions** f_m, f_h (primary)  
2. **Replace the stability functions** f_m, f_h (primary)  
3. **Use a Richardson-dependent turbulent Prandtl number** Pr_t(Ri_g) (secondary)  
4. **Use a Richardson-dependent turbulent Prandtl number** Pr_t(Ri_g) (secondary)  
5. **Use a geometric wrapper only as an optional compatibility layer** (tertiary)  
6. **Use a geometric wrapper only as an optional compatibility layer** (tertiary)  
Here’s why.  
**1. Replace f_m and f_h first**  
Your recent derivation already establishes a direct mapping  
\zeta=\frac{Ri_g}{1-5Ri_g},  
which implies  
f_m(Ri_g)=\left(1-5Ri_g\right)^2,  
for the linear Businger-Dyer form with \beta_m=5.  
This is much stronger than simply proposing another empirical curve.  
It means the Richardson-number closure is **analytically implied** by MOST itself after eliminating \zeta.  
It means the Richardson-number closure is **analytically implied** by MOST itself after eliminating \zeta.  
Instead of treating  
f(Ri)=\frac1{1+bRi}  
or  
e^{-cRi}  
as curve fits, you now have a closure derived from the similarity equations.  
That naturally suggests allowing  
f_m(Ri)\neq f_h(Ri),  
rather than forcing identical decay.  
   
⸻  
   
**2. Variable Pr_t is the natural extension**  
Once momentum and heat are separated,  
K_m=l^2S\,f_m(Ri),  
K_h=l^2S\,f_h(Ri),  
then  
Pr_t=\frac{K_m}{K_h} =\frac{f_m}{f_h}.  
Notice something important.  
Instead of **prescribing**  
Pr_t=0.85+2Ri,  
or  
Pr_t=0.85e^{Ri/0.2}+0.15,  
you can derive  
Pr_t(Ri) = \frac{f_m(Ri)}{f_h(Ri)}.  
This is considerably cleaner.  
The turbulent Prandtl number becomes a consequence of the closure rather than another empirical tuning function.  
That fits your overall philosophy of minimizing empirical parameters.  
   
⸻  
   
**3. The manifold wrapper is a higher-level correction**  
The multiplicative correction  
K_{\rm corr} = K_{\rm MOST}\, \mathcal C(d_{\rm fold})  
is attractive because it leaves existing codes untouched.  
However, from the GSPT perspective it has a different interpretation.  
It should not replace the stability function.  
Instead,  
\mathcal C  
acts as a **topological correction** accounting for proximity to a loss of normal hyperbolicity.  
In other words,  
K = l^2S f(Ri) \underbrace{\mathcal C(d_{\rm fold})}_{\text{slow-manifold correction}}.  
The Richardson closure represents **local turbulence physics**.  
The manifold correction represents **global state-space geometry**.  
The manifold correction represents **global state-space geometry**.  
Those are different pieces of the model.  
   
⸻  
   
## I would slightly modify the Julia implementation  
Rather than  
```
fm = exp(-2.5*Rig)

Pr_t = Pr0 + 2Rig

Kh = Km/Pr_t

```
I’d separate the closures completely.  
```
fm = momentum_stability(Rig)
fh = heat_stability(Rig)

Km = l_mix^2 * shear * fm
Kh = l_mix^2 * shear * fh

Pr_t = Km / Kh

```
Now  
* momentum_stability  
* heat_stability  
can be swapped independently.  
MOST, SHEBA, Grachev, your new analytical mapping, or future GSPT closures all use the same interface.  
   
⸻  
   
## Architecturally  
Within GeoABL, I’d structure the closure as  
```
Mixing Length
      │
      ▼
Ri_g
      │
      ▼
Momentum Stability Function
      │
      ▼
Heat Stability Function
      │
      ▼
Km, Kh
      │
      ▼
(Optional)
Geometric Suppression
      │
      ▼
Final Eddy Diffusivities

```
This keeps:  
* **Level 1:** local similarity physics (MOST-derived closure),  
* **Level 2:** transport anisotropy through Pr_t,  
* **Level 3:** GSPT/topological corrections,  
cleanly separated.  
That organization also matches the modular architecture you’ve been building for GeoABL, where local turbulence closure and slow-manifold geometry remain distinct components rather than being blended into a single empirical function.  
