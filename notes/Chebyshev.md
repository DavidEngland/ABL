To transform the Businger-Dyer (BD) profile for heat into a Gegenbauer form using the change of variables $t = \sqrt{b_h \zeta}$, we essentially shift from a power series in $\zeta$ to an expansion in a coordinate that linearizes the internal structure of the stability function.
1. The Target Expression
The unstable BD profile for heat is:
$$\phi_h(\zeta) = (1 - b_h \zeta)^{-1/2}$$
Substituting $t = \sqrt{b_h \zeta}$ (or $b_h \zeta = t^2$), the expression becomes:
$$\phi_h(t) = (1 - t^2)^{-1/2}$$
2. Identifying the Gegenbauer Generating Function
The generating function for Gegenbauer polynomials $C_n^{(\lambda)}(t)$ is given by:
$$(1 - 2tx + x^2)^{-\lambda} = \sum_{n=0}^{\infty} C_n^{(\lambda)}(t) x^n$$
To match our profile $(1 - t^2)^{-1/2}$ to this form, we look for a specific case. Note that if we set $x = t$ in the generating function, the expression becomes $(1 - 2t^2 + t^2)^{-\lambda} = (1 - t^2)^{-\lambda}$.
For the heat profile, we have an exponent of $-1/2$, so we set $\lambda = 1/2$.
3. The Transformation to Gegenbauer Form
When $\lambda = 1/2$, the Gegenbauer polynomials $C_n^{(1/2)}(t)$ are exactly the Legendre Polynomials$P_n(t)$. The generating function becomes:
$$(1 - 2tx + x^2)^{-1/2} = \sum_{n=0}^{\infty} P_n(t) x^n$$
By setting $x = t$ in this expansion:
$$\phi_h(t) = (1 - t^2)^{-1/2} = \sum_{n=0}^{\infty} P_n(t) t^n$$
This gives you the "Gegenbauer form" (specifically the Legendre case) of the heat profile:
$$\phi_h(\sqrt{b_h \zeta}) = \sum_{n=0}^{\infty} P_n(\sqrt{b_h \zeta}) \cdot (\sqrt{b_h \zeta})^n$$

4. Why this specific $t = \sqrt{b_h \zeta}$ mapping?
1. Singularity Alignment: The BD profile has a singularity at $b_h \zeta = 1$. By using $t = \sqrt{b_h \zeta}$, the domain of the physical stability space is mapped onto the standard interval of orthogonal polynomials ($t \in [0, 1]$).
2. Symmetry: The form $(1-t^2)^{-1/2}$ is the weight function for Chebyshev polynomials of the first kind($T_n$). This implies that the heat profile is actually the "inverse" of the very basis functions often used to solve the momentum equations.
3. Numerical Efficiency: In stable boundary layer modeling, evaluating the profile in this form allows you to use the Gegenbauer recurrence relation mentioned earlier. Instead of computing square roots and divisions repeatedly, you can propagate the profile values using simple additions and multiplications of $P_n(t)$.
5. Contrast with Momentum
If you were to do this for the momentum profile $\phi_m = (1 - b_m \zeta)^{-1/4}$, you would follow the same substitution $t = \sqrt{b_m \zeta}$, but your Gegenbauer parameter would be $\lambda = 1/4$:
$$\phi_m(t) = (1 - t^2)^{-1/4} = \sum_{n=0}^{\infty} C_n^{(1/4)}(t) t^n$$
This confirms why the momentum operator $\mathcal{L}_m$ requires that specific $1/4$ power singular weight—it is the natural orthogonal basis for the momentum profile's geometric structure.

There are some significant mathematical errors in this derivation that need to be corrected before the framework can be trusted.
The core error: the generating function substitution
The Gegenbauer generating function is:
￼
Setting ￼ gives:
￼
So far so good — but this is only valid when ￼ and ￼ simultaneously, meaning ￼. More critically, the expansion

\sum_{n=0}^\infty P_n(t), t^n

is not a standard orthogonal polynomial expansion — it’s a power series with polynomial coefficients that are themselves functions of the expansion variable. This is not what orthogonal polynomial expansions mean. A proper expansion in Legendre polynomials would be:

f(t) = \sum_{n=0}^\infty a_n P_n(t), \quad a_n = \frac{2n+1}{2}\int_{-1}^1 f(t) P_n(t), dt

with constant coefficients ￼, not ￼.
What the substitution actually gives you
After ￼, you have ￼, which is indeed the Chebyshev weight kernel — but the correct expansion is in Chebyshev polynomials of the first kind:

(1-t^2)^{-1/2} = \frac{\pi}{2}\left[1 + \frac{2}{\pi}\sum_{n=1}^\infty T_{2n}(t) \cdot \frac{(-1)^{n+1}(2n-1)!!}{(2n)!!}\cdots\right]

or more cleanly via the known generating function

\frac{1}{\sqrt{1-t^2}} = \sum_{n=0}^\infty \varepsilon_n^{-1} T_n(t)^2 \quad \text{(Neumann expansion)}

where ￼, ￼ for ￼.
The momentum case has the same problem
￼ does not equal ￼ in any orthogonally meaningful sense. The ￼ trick just recovers the generating function evaluated at a specific point — it gives you a numerical identity, not an expansion in an orthogonal basis.
What is salvageable
The physical intuition is sound:
 • The substitution ￼ does map the singularity of ￼ to ￼, aligning with the natural domain ￼ of classical orthogonal polynomials
 • ￼ genuinely is the Chebyshev-￼ weight function, so Chebyshev expansion is the natural basis for the heat profile
 • ￼ is the natural weight for Gegenbauer polynomials with ￼, making those the proper basis for momentum
The correct route is to compute the actual projection coefficients ￼ and expand in those bases — which gives you the numerical efficiency and recurrence relations you’re after, but through legitimate orthogonal projection rather than the generating-function shortcut.​​​​​​​​​​​​​​​​
