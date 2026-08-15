We consider

\phi_h(\zeta)=\frac{1}{\sqrt{1-16\zeta}},


a prototype heat–shear similarity function whose analytic structure is governed by a single branch point at \(\zeta=1/16\). The representation naturally yields a central‑binomial expansion in the neutral regime and a clean unstable‑boundary‑layer (UBL) asymptotic describing the approach to strong convection.

Neutral-regime expansion

For \(|\zeta|<1/16\), the classical binomial identity

(1-x)^{-1/2}=\sum_{n=0}^{\infty}\frac{\binom{2n}{n}}{4^n}x^n,\qquad x=16\zeta,


gives the exact series

\phi_h(\zeta)=\sum_{n=0}^{\infty}\binom{2n}{n}4^n\zeta^n
=1+8\zeta+96\zeta^2+1280\zeta^3+17920\zeta^4+\cdots.


This central‑binomial expansion provides closed‑form coefficients, a simple recurrence, and a transparent radius of convergence set by the physical singularity at \(\zeta=1/16\). It is the natural analytic basis for neutral‑regime calibration and for constructing rational or piecewise extensions.

Unstable (UBL) asymptotic: strong‑convective limit

For \(\zeta<0\), write \(\eta=-\zeta>0\). Then

\phi_h(\zeta)=\frac{1}{\sqrt{1+16\eta}}
\sim \frac{1}{4\sqrt{\eta}}
\left(1-\frac{1}{32\eta}+\frac{3}{2048\eta^2}-\cdots\right),
\qquad \eta\to\infty.


This is the UBL heat–shear asymptotic, and it is the principal analytic result:

• the strong‑convective limit is algebraic, \(\phi_h\propto |\zeta|^{-1/2}\);
• corrections form a regular inverse‑power hierarchy;
• no oscillatory or exponential behavior appears.


This asymptotic form provides a physically interpretable scaling for strongly unstable boundary layers and a natural anchor for matched‑regime closures.

Stable and very stable regimes

For \(0<\zeta<1/16\), \(\phi_h\) increases monotonically and diverges at the branch point. Beyond \(\zeta=1/16\) the expression becomes complex. Thus the raw form cannot serve as a universal stable‑regime closure; practical implementations require a regularized, capped, or rational continuation that preserves continuity with the neutral expansion.

Numerical considerations

Evaluation via

\phi_h=\exp\!\left(-\frac12\log(1-16\zeta)\right)


with log1p(-16ζ) avoids cancellation near neutrality. The coefficients

c_{n+1}=c_n\,\frac{8(2n+1)}{n+1},\qquad c_0=1,


enable efficient generation of truncated neutral‑regime series, and Horner evaluation minimizes roundoff.

Summary

The function \(\phi_h(\zeta)=(1-16\zeta)^{-1/2}\) admits a central‑binomial neutral expansion and a UBL asymptotic

\phi_h(\zeta)\sim \frac{1}{4\sqrt{-\zeta}},\qquad \zeta\to -\infty,


which together provide a complete analytic description of the unstable and near‑neutral regimes. The positive‑\(\zeta\) singularity at \(1/16\) limits direct applicability in stable conditions, motivating rational or piecewise continuations in practical boundary‑layer closures.

---

