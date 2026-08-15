# 🔍 Why `arcoth`, `artanh`, `ln`, `Log`, and Complex Transformations Matter

## 1. Unified Structure: Hyperbolic Inverses and Logarithms

```math
\text{arcoth}(x) = \frac{1}{2} \ln\left(\frac{x + 1}{x - 1}\right)
\quad
\text{artanh}(x) = \frac{1}{2} \ln\left(\frac{1 + x}{1 - x}\right)
```

- These are logarithmic in nature.
- `ln(x)` and `Log(z)` serve as the analytic backbone.

---

## 2. Complex Transformations: Link to Trigonometric Functions

```math
\text{artanh}(i x) = i \cdot \arctan(x)
\quad
\text{arcoth}(i x) = i \cdot \text{arccot}(x)
```

- Hyperbolic ↔ Trigonometric via complex substitution.
- Logarithmic definitions extend naturally to the complex plane.

---

## 3. Series Representations

```math
\text{artanh}(x) = \sum_{n=0}^{\infty} \frac{x^{2n+1}}{2n+1}, \quad |x| < 1
```

```math
\text{arcoth}(x) = \sum_{n=0}^{\infty} \frac{x^{-(2n+1)}}{2n+1}, \quad |x| > 1
```

- These mirror `arctan(x)` and `arccot(x)` under `x → i·x`.

---

## 4. Connection to Exponential and Bernoulli Numbers

```math
\sinh(x) = \frac{e^x - e^{-x}}{2}
\quad
\cosh(x) = \frac{e^x + e^{-x}}{2}
```

Bernoulli expansions:

```math
x \cdot \cot(x) = 1 - 2 \sum_{n=1}^{\infty} \frac{(-1)^n B_{2n}}{(2n)!} (2x)^{2n}
```

```math
x \cdot \coth(x) = 1 + 2 \sum_{n=1}^{\infty} \frac{B_{2n}}{(2n)!} (2x)^{2n}
```

- Lagrange inversion yields series for `arccot(x)` and `arcoth(x)`.

---

## 5. Analytic Continuation and Branch Cuts

- `Log(z)` enables analytic continuation.
- Branch cuts governed by logarithmic and square-root structure.
- Crucial in complex analysis, physics, and number theory.

---

## 6. Summary Table

| Function     | Logarithmic Form                          | Series Expansion                     | Complex/Trig Link                     | Bernoulli Connection |
|--------------|--------------------------------------------|--------------------------------------|----------------------------------------|----------------------|
| `artanh(x)`  | `(1/2)·ln((1 + x)/(1 - x))`               | `Σ x^(2n+1)/(2n+1)`                  | `artanh(i·x) = i·arctan(x)`            | ✓                    |
| `arcoth(x)`  | `(1/2)·ln((x + 1)/(x - 1))`               | `Σ x^-(2n+1)/(2n+1)`                 | `arcoth(i·x) = i·arccot(x)`            | ✓                    |
| `ln(x)`      | —                                          | —                                    | —                                      | —                    |
| `cot(x)`     | —                                          | —                                    | —                                      | Bernoulli expansion  |
| `coth(x)`    | —                                          | —                                    | —                                      | Bernoulli expansion  |

---

## 7. Why This Is a Big Deal

- Unifies hyperbolic, trigonometric, logarithmic, and exponential functions.
- Enables series expansions, special values, and analytic continuation.
- Bernoulli numbers and zeta values encode deep arithmetic structure.
- Applications in analysis, number theory, physics, and engineering.

---

## 8. The Logarithm via Inverse Hyperbolic Functions

```math
\ln(x) = \text{artanh}\left(\frac{x^2 - 1}{x^2 + 1}\right)
       = \text{arsinh}\left(\frac{x^2 - 1}{2x}\right)
       = \text{arcoth}\left(\frac{x^2 + 1}{x^2 - 1}\right)
```

- Each real or complex `x` corresponds to a unique hyperbolic angle.
- Use `artanh(x)` for `|x| < 1`, `arcoth(x)` for `|x| > 1`.

Reciprocal identity:

```math
\text{artanh}(x) = \text{arcoth}\left(\frac{1}{x}\right), \quad x \ne 0
```

---

## 📚 References

- *Abramowitz & Stegun*, Handbook of Mathematical Functions
- *NIST DLMF*, Digital Library of Mathematical Functions
- Wikipedia: Inverse hyperbolic functions, Bernoulli numbers, Lagrange inversion theorem

---
