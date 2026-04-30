Subject: Fractional dimension interpretation for momentum and scalar MOST exponents

Dick and Arastoo,

I wanted to float a simple way to think about the fractional-dimension idea in the ultraspherical/Gegenbauer framing.

The key point is that there are really two parameters in play, and I think we should keep them notationally separate:

1. The Gegenbauer parameter, which I will call $\Lambda_q$.
2. The profile exponent parameter in a closure of the form
$$
\phi_q(\zeta)=a_q(1-b_q\zeta)^{-1/\lambda_q}.
$$

These are related by
$$
\Lambda_q = \frac{1}{\lambda_q}.
$$

For the ultraspherical family, the effective dimension is tied to the Gegenbauer parameter through
$$
\Lambda_q = \frac{d_q-2}{2},
$$
so equivalently
$$
d_q = 2\Lambda_q + 2 = 2 + \frac{2}{\lambda_q}.
$$

This gives a very simple interpretation of the classical unstable MOST exponents:

- Momentum: exponent $-1/4$ means $\Lambda_m = 1/4$, hence
$$
d_m = 2 + 2\left(\frac{1}{4}\right)=2.5.
$$
- Heat: exponent $-1/2$ means $\Lambda_h = 1/2$, hence
$$
d_h = 2 + 2\left(\frac{1}{2}\right)=3.
$$

So one possible interpretation is:

- the heat/scalar branch looks fully three-dimensional,
- while momentum behaves as if the turbulence is geometrically constrained to an effective dimension of about $2.5$ in the surface layer.

That feels at least qualitatively consistent with momentum being more strongly constrained by wall-normal suppression and anisotropy than scalar transport.

The part I am less certain about, and where I would value your thoughts, is how far to push this for other tracers.

For humidity, CO2, methane, aerosols, etc., I do not think we can yet claim universal canonical exponents analogous to the momentum and heat values. My current view is:

1. Conservative tracers near neutral conditions probably start near the scalar limit $d_q \approx 3$.
2. Nonconservative tracers may depart from that because of chemistry, phase change, deposition, canopy exchange, or particle settling.
3. In that sense, $d_q$ may be better treated as a tracer-dependent inferred quantity rather than a universal constant.

This suggests a possible next step:

1. Keep the classical momentum and heat values as anchor points.
2. Write $d_q$ (or equivalently $\Lambda_q$) as a dynamic or fitted quantity for additional tracers.
3. See whether observed tracer transfer functions collapse better when plotted against an inferred $d_q$.

If this line of thinking seems reasonable, I can write up a short note with the derivation and show how it connects to anisotropy arguments and to the usual MOST forms.

David
