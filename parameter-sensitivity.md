We’ll expand the parameter-sensitivity derivatives
\frac{\partial f}{\partial \alpha_h}, \quad \frac{\partial f}{\partial \beta_m}, \quad \text{etc.}
for a standard Hodgkin–Huxley–type gating formulation, where each gate variable x \in \{m, h, n, \dots\} satisfies:

\frac{dx}{dt} = \alpha_x(V) (1 - x) - \beta_x(V) x.

⸻

1. General setup

Let
f_x(V,x;\alpha_x,\beta_x) = \alpha_x(V)(1 - x) - \beta_x(V)x
be the right-hand side of the gating ODE.

Then sensitivities to the kinetic parameters are:

\frac{\partial f_x}{\partial \alpha_x} = (1 - x),
\qquad
\frac{\partial f_x}{\partial \beta_x} = -x.

If \alpha_x or \beta_x themselves depend on parameters (e.g. rate constants or slope factors), then we expand further.

⸻

2. Typical parameterized form of rate functions

For the original HH sodium activation/inactivation and potassium activation gates:

\alpha_m(V) = \frac{A_m (V + B_m)}{1 - e^{-(V + B_m)/C_m}}, \quad
\beta_m(V) = D_m e^{-(V + E_m)/F_m},

\alpha_h(V) = G_h e^{-(V + H_h)/I_h}, \quad
\beta_h(V) = \frac{J_h}{1 + e^{-(V + K_h)/L_h}},

(and similar for n).

Then each parameter-sensitivity derivative like \partial f_x / \partial A_m expands by the chain rule:

\frac{\partial f_x}{\partial A_m}
= \frac{\partial f_x}{\partial \alpha_x} \frac{\partial \alpha_x}{\partial A_m}
	•	\frac{\partial f_x}{\partial \beta_x} \frac{\partial \beta_x}{\partial A_m}.

Since \beta_x is independent of A_m, the second term is zero, giving

\frac{\partial f_m}{\partial A_m}
= (1 - m)\, \frac{\partial \alpha_m}{\partial A_m}.

⸻

3. Explicit symbolic derivatives for typical rates

Let’s compute them symbolically.

(a) For \alpha_m(V) = \dfrac{A_m (V + B_m)}{1 - e^{-(V + B_m)/C_m}}

\begin{aligned}
\frac{\partial \alpha_m}{\partial A_m}
&= \frac{V + B_m}{1 - e^{-(V + B_m)/C_m}}, \\[6pt]
\frac{\partial \alpha_m}{\partial B_m}
&= \frac{A_m}{1 - e^{-(V + B_m)/C_m}}
	•	A_m (V + B_m)\,
\frac{e^{-(V + B_m)/C_m}}{C_m(1 - e^{-(V + B_m)/C_m})^2}, \\[6pt]
\frac{\partial \alpha_m}{\partial C_m}
&= A_m (V + B_m)\,
\frac{e^{-(V + B_m)/C_m}(V + B_m)/C_m^2}
{(1 - e^{-(V + B_m)/C_m})^2}.
\end{aligned}

Then the propagated sensitivities are:

\boxed{
\begin{aligned}
\frac{\partial f_m}{\partial A_m} &= (1 - m)\,\frac{V + B_m}{1 - e^{-(V + B_m)/C_m}}, \\[4pt]
\frac{\partial f_m}{\partial B_m} &= (1 - m)\left[\frac{A_m}{1 - e^{-(V + B_m)/C_m}} +
A_m (V + B_m)\frac{e^{-(V + B_m)/C_m}}{C_m(1 - e^{-(V + B_m)/C_m})^2}\right], \\[4pt]
\frac{\partial f_m}{\partial C_m} &= (1 - m)A_m (V + B_m)
\frac{e^{-(V + B_m)/C_m}(V + B_m)/C_m^2}
{(1 - e^{-(V + B_m)/C_m})^2}.
\end{aligned}
}

⸻

(b) For \beta_m(V) = D_m e^{-(V + E_m)/F_m}

\frac{\partial \beta_m}{\partial D_m} = e^{-(V + E_m)/F_m}, \quad
\frac{\partial \beta_m}{\partial E_m} = -\frac{D_m}{F_m} e^{-(V + E_m)/F_m}, \quad
\frac{\partial \beta_m}{\partial F_m} = D_m e^{-(V + E_m)/F_m} \frac{V + E_m}{F_m^2}.

Propagating to f_m:

\boxed{
\begin{aligned}
\frac{\partial f_m}{\partial D_m} &= -m\, e^{-(V + E_m)/F_m}, \\[4pt]
\frac{\partial f_m}{\partial E_m} &= m\,\frac{D_m}{F_m}\, e^{-(V + E_m)/F_m}, \\[4pt]
\frac{\partial f_m}{\partial F_m} &= -m\,D_m\, e^{-(V + E_m)/F_m}\, \frac{V + E_m}{F_m^2}.
\end{aligned}
}

⸻

(c) For \alpha_h(V) = G_h e^{-(V + H_h)/I_h}

\boxed{
\begin{aligned}
\frac{\partial f_h}{\partial G_h} &= (1 - h) e^{-(V + H_h)/I_h}, \\[4pt]
\frac{\partial f_h}{\partial H_h} &= (1 - h)\,G_h\,\frac{e^{-(V + H_h)/I_h}}{I_h}, \\[4pt]
\frac{\partial f_h}{\partial I_h} &= -(1 - h)\,G_h\, e^{-(V + H_h)/I_h}\frac{V + H_h}{I_h^2}.
\end{aligned}
}

⸻

(d) For \beta_h(V) = \dfrac{J_h}{1 + e^{-(V + K_h)/L_h}}

\begin{aligned}
\frac{\partial \beta_h}{\partial J_h} &= \frac{1}{1 + e^{-(V + K_h)/L_h}}, \\[4pt]
\frac{\partial \beta_h}{\partial K_h} &= J_h \frac{e^{-(V + K_h)/L_h}}{L_h(1 + e^{-(V + K_h)/L_h})^2}, \\[4pt]
\frac{\partial \beta_h}{\partial L_h} &= J_h (V + K_h)
\frac{e^{-(V + K_h)/L_h}}{L_h^2(1 + e^{-(V + K_h)/L_h})^2}.
\end{aligned}

Then:

\boxed{
\begin{aligned}
\frac{\partial f_h}{\partial J_h} &= -h\,\frac{1}{1 + e^{-(V + K_h)/L_h}}, \\[4pt]
\frac{\partial f_h}{\partial K_h} &= -h\, J_h \frac{e^{-(V + K_h)/L_h}}{L_h(1 + e^{-(V + K_h)/L_h})^2}, \\[4pt]
\frac{\partial f_h}{\partial L_h} &= -h\, J_h (V + K_h)
\frac{e^{-(V + K_h)/L_h}}{L_h^2(1 + e^{-(V + K_h)/L_h})^2}.
\end{aligned}
}

⸻

4. Compact matrix form (optional)

For each gate x:

\frac{\partial f_x}{\partial \mathbf{p}_x} =
\begin{bmatrix}
(1 - x)\frac{\partial \alpha_x}{\partial \mathbf{p}_x} -
x\frac{\partial \beta_x}{\partial \mathbf{p}_x}
\end{bmatrix},

where \mathbf{p}_x = (A_x, B_x, C_x, D_x, E_x, F_x, \dots).

⸻

Would you like me to extend this to second-order sensitivities
(e.g. \partial^2 f_x / \partial \alpha_x \partial V, useful for adjoint or Fisher information analysis)?