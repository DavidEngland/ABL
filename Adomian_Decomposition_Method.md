# Adomian Decomposition Method (ADM) — Practical Reference

Purpose  
Provide a concise, implementation-oriented explanation of the Adomian Decomposition Method (ADM) for solving nonlinear differential equations and integral equations. Intended for readers with strong mathematics/physics background who want a ready-to-run recipe and practical tips (including applicability to boundary-layer problems).

Prerequisites
- Ordinary/partial differential equations, linear operators
- Familiarity with series expansions and operator inverses (Green's functions or integration operators)
- Basic symbolic algebra (helpful for generating Adomian polynomials)

1. Overview and when to use ADM
- ADM constructs a solution as a rapidly convergent series u = Σ_{n=0}^∞ u_n without needing linearization, perturbation parameters, or small/large assumptions.
- Useful for nonlinear ODEs/PDEs, integral equations, and cases where classical perturbation fails or a closed-form inverse of a dominant linear operator exists.
- Cautions: convergence must be checked; for stiff or strongly singular problems additional transforms or accelerations (Padé, Wynn-epsilon) may be required.

2. Operator formulation (core idea)
Given a nonlinear problem in operator form
L[u] + R[u] + N[u] = g,                        (1)
where
- L is the highest-order linear operator (invertible; choose it),
- R is the remainder linear operator (lower order),
- N is the nonlinear operator,
- g is the source / inhomogeneity / boundary data contribution.

Apply the inverse operator L^{-1} (integration or Green's operator) to obtain
u = L^{-1}[g] − L^{-1}[R[u]] − L^{-1}[N[u]] + terms from boundary conditions.    (2)

3. Decomposition ansatz
Decompose the solution and nonlinearity as
u = Σ_{k=0}^∞ u_k,                               (3)
N[u] = Σ_{k=0}^∞ A_k,                             (4)
where the A_k are the Adomian polynomials dependent on {u_0, ..., u_k}.

4. Construction of Adomian polynomials
General generating formula:
A_k = (1/k!) (d^k/dλ^k) [ N( Σ_{m=0}^∞ λ^m u_m ) ]_{λ=0}.

Practical rules for common nonlinearities:
- For N[u] = u^2:
  A_0 = u_0^2,
  A_1 = 2 u_0 u_1,
  A_2 = 2 u_0 u_2 + u_1^2,
  A_3 = 2 u_0 u_3 + 2 u_1 u_2, ...
- For N[u] = exp(u): generate derivatives of exp(Σ λ^m u_m) and evaluate at λ=0.
- For composite or complicated N, symbolic differentiation (sympy/Maple/Mathematica) is easiest.

5. Recursive scheme
Choose u_0 to satisfy the linear inhomogeneous part and boundary/initial conditions:
u_0 = L^{-1}[g] + (terms matching boundary/initial conditions).

Then iterate for n ≥ 0:
u_{n+1} = − L^{-1}[ R[u_n] + A_n ].                 (5)

Thus each u_{n+1} is obtained by applying L^{-1} to known terms (previous u_j and A_j). Sum partial series U^{(N)} = Σ_{k=0}^{N} u_k is the N-term ADM approximation.

6. Worked ODE example (illustrative, fully worked to low order)
Solve simple initial-value problem (toy): u'(x) = u(x)^2 + x, u(0)=0.
- Write as L[u] = u', so L^{-1}[·] = ∫_0^x (·) ds and boundary condition u(0)=0 imposes integration constant.
- Decompose u = Σ u_k, N[u] = u^2 = Σ A_k.

Start:
u_0(x) = L^{-1}[x] = ∫_0^x s ds = x^2/2.

Compute A_0 = u_0^2 = (x^2/2)^2 = x^4/4.

Then
u_1 = − L^{-1}[ A_0 ] = − ∫_0^x (s^4/4) ds = − s^5/(20) |_0^x = − x^5/20.

Next A_1 = 2 u_0 u_1 = 2*(x^2/2)*(−x^5/20) = − x^7/20.

u_2 = − L^{-1}[ A_1 ] = + ∫_0^x s^7/20 ds = + x^8/(160).

So partial sum up to second correction:
u ≈ u_0 + u_1 + u_2 = x^2/2 − x^5/20 + x^8/160 + O(x^{11}).

Remarks:
- This matches the Maclaurin expansion of the true solution near x=0.
- The procedure is algebraic and can be continued to desired order.

7. Boundary-layer / PDE remarks
- For boundary-layer problems (e.g., Blasius-like or nonlinear diffusion), choose L to be the dominant linear differential operator (e.g., highest derivative ∂^n/∂x^n) and ensure L^{-1} respects boundary conditions (Green's function or integral operator).
- If domain is semi-infinite (0,∞) with conditions at infinity, ADM series may require re-summation or asymptotic matching to respect far-field BC; common remedies:
  - use variable transforms (compactify domain),
  - introduce composite expansions (ADM for inner region + matched outer solution),
  - use Padé approximants for series acceleration.

8. Convergence, error control, and acceleration
- Convergence is problem-dependent. Practical checks:
  - monitor residual R_N(x) = L[u^{(N)}] + R[u^{(N)}] + N[u^{(N)}] − g; stop when ||R_N|| below tolerance.
  - inspect term norms ||u_{n+1}||/||u_n|| for geometric decay.
- Acceleration techniques:
  - Padé approximants of partial sums often extend radius of validity.
  - Wynn-epsilon or Shanks transforms can accelerate slowly convergent series.
  - Homotopy-Padé and parameter embeddings combine ADM with continuation to improve convergence.

9. Implementation tips
- Generate Adomian polynomials symbolically for nontrivial N using a CAS; store recurrence formulas.
- If L^{-1} is an integral operator, precompute integrals of basis monomials to speed up term evaluation.
- Use automatic differentiation or symbolic algebra to compute A_k when N is complicated (e.g., functions of derivatives).
- For PDEs discretize transverse coordinates if needed and apply ADM in the primary coordinate (method of lines + ADM in the stiff dimension).

10. Practical checklist before using ADM
- Ensure an explicit, bounded inverse L^{-1} exists and implementable with BCs.
- Check local analyticity of N[u] in the function-space region of interest (Adomian generation requires differentiability).
- If the domain includes infinity, plan transforms or acceleration.
- Validate against a numerical integrator (Runge–Kutta / spectral method) for at least one case.

11. References and further reading (suggested keywords)
- "Adomian Decomposition Method" (original papers and review articles)
- "Adomian polynomials" + "Padé acceleration"
- Use CAS implementations (SymPy examples: `adomian` recipe) for practice.

Appendix: quick pseudocode (conceptual)
```
Given L, R, N, g and BCs:
1. compute u0 = L^{-1}[g] + BC_terms
2. for n = 0 .. N-1:
     compute A_n from u_0..u_n
     compute rhs = R[u_n] + A_n
     u_{n+1} = − L^{-1}[rhs]
3. return partial_sum = Σ_{k=0}^N u_k
```

Notes for boundary-layer community
- ADM produces analytic-like series that can illuminate asymptotic structure of near-wall layers.
- For strongly nonlinear SBL models with history dependence or turbulence closures, ADM can produce semi-analytic kernels useful for reduced-order modeling, but caution is needed for long-range (in z) behavior and for ensuring physically consistent far-field BCs.

Contact / usage
- Implement small examples in symbolic Python/Mathematica to build intuition.
- Use ADM for quick analytic approximations, initial guesses for shooting/continuation methods, and for testing asymptotic balances before heavy numerics.
